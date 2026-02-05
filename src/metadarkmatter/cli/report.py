"""
Report command for generating unified HTML reports.

Creates self-contained HTML reports combining all visualization types
with interactive features and professional styling.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from metadarkmatter.cli.utils import QuietConsole, extract_sample_name, spinner_progress
from metadarkmatter.core.io_utils import read_dataframe
from metadarkmatter.core.metadata import GenomeMetadata

app = typer.Typer(
    name="report",
    help="Generate unified HTML reports from classification results",
    no_args_is_help=True,
)

console = Console()


@app.command(name="generate")
def generate_report(
    classifications: Path = typer.Option(
        ...,
        "--classifications", "-c",
        help="Classification results file (CSV or Parquet)",
        exists=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output HTML report path",
    ),
    sample_name: str = typer.Option(
        None,
        "--sample-name", "-n",
        help="Sample name for report header (defaults to input filename)",
    ),
    title: str = typer.Option(
        "Metadarkmatter Classification Report",
        "--title", "-t",
        help="Report title",
    ),
    ani_matrix: Path | None = typer.Option(
        None,
        "--ani", "-a",
        help="ANI matrix file (CSV) for heatmap visualization",
        exists=True,
        dir_okay=False,
    ),
    aai_matrix: Path | None = typer.Option(
        None,
        "--aai",
        help="AAI matrix file (CSV) for genus-level heatmap visualization",
        exists=True,
        dir_okay=False,
    ),
    metadata: Path | None = typer.Option(
        None,
        "--metadata", "-m",
        help="Genome metadata file (TSV) for species-level breakdown",
        exists=True,
        dir_okay=False,
    ),
    recruitment_data: Path | None = typer.Option(
        None,
        "--recruitment", "-r",
        help="Recruitment data file (CSV) for recruitment plots",
        exists=True,
        dir_okay=False,
    ),
    bam: Path | None = typer.Option(
        None,
        "--bam", "-b",
        help="BAM file for generating recruitment data on-the-fly",
        exists=True,
        dir_okay=False,
    ),
    theme: str = typer.Option(
        "light",
        "--theme",
        help="Color theme: 'light' or 'dark'",
    ),
    alignment_mode: str = typer.Option(
        "nucleotide",
        "--alignment-mode",
        help="Alignment mode: 'nucleotide' (uses ANI) or 'protein' (uses AAI)",
    ),
    max_scatter_points: int = typer.Option(
        50000,
        "--max-points",
        help="Maximum points in scatter plots (for performance)",
        min=1000,
    ),
    max_table_rows: int = typer.Option(
        10000,
        "--max-table-rows",
        help="Maximum rows in data table (for file size)",
        min=100,
    ),
    max_phylo_clusters: int = typer.Option(
        20,
        "--max-phylo-clusters",
        help="Maximum novel clusters in phylogenetic context heatmap",
        min=5,
        max=100,
    ),
    max_phylo_references: int = typer.Option(
        50,
        "--max-phylo-references",
        help="Maximum reference genomes in phylogenetic context heatmap",
        min=10,
        max=200,
    ),
    tree: Path | None = typer.Option(
        None,
        "--tree",
        help="Newick tree file for phylogenetic context. "
             "If not provided, neighbor-joining tree is built from ANI matrix.",
        exists=True,
        dir_okay=False,
    ),
    no_phylogeny: bool = typer.Option(
        False,
        "--no-phylogeny",
        help="Skip phylogeny tab generation (faster for large datasets).",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Generate a unified HTML report from classification results.

    Creates a self-contained, interactive HTML report with multiple tabs:

    - Overview: Summary metrics and classification distribution charts
    - Distributions: Histograms and 2D novelty vs uncertainty scatter plot
    - Recruitment: Read recruitment plots (if BAM/data provided)
    - Species: Species-level breakdown charts (if metadata provided)
    - Genomes: Per-genome breakdown charts
    - ANI Matrix: Heatmap of genome-genome ANI values (if provided)
    - Phylogeny: Interactive phylogenetic tree with novel cluster placement (if ANI provided)
    - Data: Interactive searchable/sortable classification table

    Example:

        # Basic report from classification results
        metadarkmatter report generate \\
            --classifications results.csv \\
            --output report.html \\
            --sample-name "Sample_001"

        # Report with species-level breakdown
        metadarkmatter report generate \\
            --classifications results.csv \\
            --metadata genome_metadata.tsv \\
            --output report.html \\
            --sample-name "Sample_001"

        # Full report with ANI matrix, recruitment plots, and species
        metadarkmatter report generate \\
            --classifications results.csv \\
            --output report.html \\
            --ani ani_matrix.csv \\
            --metadata genome_metadata.tsv \\
            --bam mapped.bam \\
            --sample-name "Environmental Sample"

        # Report with custom phylogenetic tree
        metadarkmatter report generate \\
            --classifications results.csv \\
            --ani ani_matrix.csv \\
            --tree custom_tree.nwk \\
            --output report.html

        # Dark theme report
        metadarkmatter report generate \\
            --classifications results.csv \\
            --output report.html \\
            --theme dark
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Report Generator[/bold blue]\n")

    # Import report modules
    try:
        from metadarkmatter.visualization.plots.base import PlotConfig, ThresholdConfig
        from metadarkmatter.visualization.report import (
            ReportConfig,
            ReportGenerator,
        )
    except ImportError as e:
        console.print(
            f"[red]Error: Required modules not available: {e}[/red]\n"
            "[dim]Ensure plotly is installed: pip install plotly[/dim]"
        )
        raise typer.Exit(code=1) from None

    # Validate theme
    if theme not in ("light", "dark"):
        console.print(
            f"[red]Error: Invalid theme '{theme}'. Use 'light' or 'dark'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Determine sample name
    if sample_name is None:
        sample_name = extract_sample_name(classifications)

    # Ensure output directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Load classification data
    out.print("[bold]Step 1:[/bold] Loading classification data...")

    with spinner_progress("Reading classifications...", console, quiet):
        try:
            df = read_dataframe(classifications)
        except Exception as e:
            console.print(f"[red]Error reading classifications: {e}[/red]")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Loaded {len(df):,} classifications[/green]")

    if verbose:
        # Show classification breakdown
        counts = df.group_by("taxonomic_call").len().sort("len", descending=True)
        for row in counts.iter_rows(named=True):
            pct = row["len"] / len(df) * 100
            out.print(f"  [dim]{row['taxonomic_call']}: {row['len']:,} ({pct:.1f}%)[/dim]")

    # Step 2: Load optional data
    ani_df = None
    if ani_matrix:
        out.print("\n[bold]Step 2a:[/bold] Loading ANI matrix...")
        try:
            ani_df = pl.read_csv(ani_matrix)
            out.print(f"  [green]Loaded ANI matrix ({len(ani_df)} genomes)[/green]")
        except Exception as e:
            console.print(f"[yellow]Warning: Could not load ANI matrix: {e}[/yellow]")

    aai_df = None
    if aai_matrix:
        out.print("\n[bold]Step 2a2:[/bold] Loading AAI matrix...")
        try:
            aai_df = pl.read_csv(aai_matrix)
            out.print(f"  [green]Loaded AAI matrix ({len(aai_df)} genomes)[/green]")
        except Exception as e:
            console.print(f"[yellow]Warning: Could not load AAI matrix: {e}[/yellow]")

    genome_metadata: GenomeMetadata | None = None
    if metadata:
        out.print("\n[bold]Step 2b:[/bold] Loading genome metadata...")
        try:
            genome_metadata = GenomeMetadata.from_file(metadata)
            out.print(
                f"  [green]Loaded metadata for {genome_metadata.genome_count} genomes "
                f"({genome_metadata.species_count} species)[/green]"
            )
        except Exception as e:
            console.print(f"[yellow]Warning: Could not load metadata: {e}[/yellow]")

    recruitment_df = None
    if recruitment_data:
        out.print("\n[bold]Step 2c:[/bold] Loading recruitment data...")
        try:
            recruitment_df = pl.read_csv(recruitment_data)
            out.print(f"  [green]Loaded {len(recruitment_df):,} recruitment records[/green]")
        except Exception as e:
            console.print(f"[yellow]Warning: Could not load recruitment data: {e}[/yellow]")
    elif bam:
        out.print("\n[bold]Step 2c:[/bold] Extracting recruitment data from BAM...")

        try:
            from metadarkmatter.external import Samtools

            if not Samtools.check_available():
                console.print(
                    "[yellow]Warning: samtools not found, skipping recruitment plots[/yellow]"
                )
            else:
                from metadarkmatter.core.recruitment import load_recruitment_data

                with spinner_progress("Reading BAM...", console, quiet):
                    recruitment_df = load_recruitment_data(
                        bam_path=bam,
                        max_records=max_scatter_points * 5,
                    )

                out.print(f"  [green]Loaded {len(recruitment_df):,} alignments[/green]")
        except Exception as e:
            console.print(f"[yellow]Warning: Could not extract from BAM: {e}[/yellow]")

    # Step 3: Generate report
    out.print("\n[bold]Step 3:[/bold] Generating report...")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console if not quiet else None,
    ) as progress:
        task = progress.add_task(description="Building visualizations...", total=None)

        try:
            config = ReportConfig(
                sample_name=sample_name,
                title=title,
                theme=theme,
                alignment_mode=alignment_mode,
                max_scatter_points=max_scatter_points,
                max_table_rows=max_table_rows,
                max_phylo_clusters=max_phylo_clusters,
                max_phylo_references=max_phylo_references,
                skip_phylogeny=no_phylogeny,
                user_tree_path=tree,
                plot_config=PlotConfig(),
                thresholds=ThresholdConfig(),
            )

            generator = ReportGenerator(
                classifications=df,
                config=config,
                recruitment_data=recruitment_df,
                ani_matrix=ani_df,
                aai_matrix=aai_df,
                genome_metadata=genome_metadata,
            )

            progress.update(task, description="Writing HTML...")
            generator.generate(output)

        except Exception as e:
            console.print(f"\n[red]Error generating report: {e}[/red]")
            if verbose:
                console.print_exception()
            raise typer.Exit(code=1) from None

    # Calculate file size
    file_size_mb = output.stat().st_size / (1024 * 1024)

    out.print(f"  [green]Report saved to {output}[/green]")
    out.print(f"  [dim]File size: {file_size_mb:.1f} MB[/dim]")

    # Summary
    out.print("\n[bold green]Report generation complete![/bold green]")
    out.print("\n[bold]Report contents:[/bold]")
    out.print(f"  - Total reads: {len(df):,}")
    out.print(f"  - ANI matrix: {'Yes' if ani_df is not None else 'No'}")
    out.print(f"  - AAI matrix: {'Yes' if aai_df is not None else 'No'}")
    out.print(f"  - Recruitment plots: {'Yes' if recruitment_df is not None else 'No'}")
    if genome_metadata:
        out.print(f"  - Species breakdown: {genome_metadata.species_count} species")
    else:
        out.print("  - Species breakdown: No (use --metadata for species section)")
    out.print(f"  - Theme: {theme}")
    out.print(f"\n[dim]Open in browser: file://{output.absolute()}[/dim]\n")


@app.command(name="multi")
def multi_sample_report(
    input_dir: Path = typer.Option(
        ...,
        "--input-dir", "-i",
        help="Directory containing classification result files",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output HTML report path",
    ),
    pattern: str = typer.Option(
        "*classifications*.csv",
        "--pattern", "-p",
        help="Glob pattern to match classification files",
    ),
    title: str = typer.Option(
        "Multi-Sample Comparison Report",
        "--title", "-t",
        help="Report title",
    ),
    theme: str = typer.Option(
        "light",
        "--theme",
        help="Color theme: 'light' or 'dark'",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Generate a multi-sample comparison report.

    Creates an HTML report comparing classification results across
    multiple samples in a directory. The report includes:

    - Summary table with statistics for each sample
    - Stacked and grouped bar charts of classifications
    - Heatmap of classification proportions
    - Sample diversity scatter plot
    - Novelty distribution box plots
    - Novel diversity trend line

    Example:

        metadarkmatter report multi \\
            --input-dir ./results/ \\
            --output comparison.html \\
            --pattern "*_classifications.csv"
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Multi-Sample Comparison Report[/bold blue]\n")

    # Find classification files
    files = sorted(input_dir.glob(pattern))

    if not files:
        console.print(
            f"[red]Error: No files found matching pattern '{pattern}' in {input_dir}[/red]"
        )
        raise typer.Exit(code=1) from None

    out.print(f"[bold]Found {len(files)} sample files:[/bold]")
    for f in files[:10]:
        out.print(f"  [dim]{f.name}[/dim]")
    if len(files) > 10:
        out.print(f"  [dim]... and {len(files) - 10} more[/dim]")

    # Validate theme
    if theme not in ("light", "dark"):
        console.print(
            f"[red]Error: Invalid theme '{theme}'. Use 'light' or 'dark'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Ensure output directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    # Load samples
    out.print("\n[bold]Step 1:[/bold] Loading sample data...")

    with spinner_progress("Loading samples...", console, quiet):
        sample_data: dict[str, pl.DataFrame] = {}
        failed = []

        for f in files:
            sample_name = extract_sample_name(f)
            try:
                df = read_dataframe(f)
                sample_data[sample_name] = df
                if verbose:
                    out.print(f"  [dim]Loaded {sample_name}: {len(df):,} reads[/dim]")
            except Exception as e:
                failed.append((f.name, str(e)))

    if failed:
        for name, err in failed:
            console.print(f"[yellow]Warning: Could not load {name}: {err}[/yellow]")

    if not sample_data:
        console.print("[red]Error: No valid classification files could be loaded[/red]")
        raise typer.Exit(code=1) from None

    out.print(f"  [green]Loaded {len(sample_data)} samples[/green]")

    # Compute total reads
    total_reads = sum(len(df) for df in sample_data.values())
    out.print(f"  [dim]Total reads across all samples: {total_reads:,}[/dim]")

    # Generate report
    out.print("\n[bold]Step 2:[/bold] Generating comparison report...")

    with spinner_progress("Building visualizations...", console, quiet):
        try:
            from metadarkmatter.visualization.report.multi_generator import (
                MultiSampleConfig,
                MultiSampleReportGenerator,
            )

            config = MultiSampleConfig(title=title, theme=theme)
            generator = MultiSampleReportGenerator(sample_data, config=config)
            generator.generate(output)

        except Exception as e:
            console.print(f"\n[red]Error generating report: {e}[/red]")
            if verbose:
                console.print_exception()
            raise typer.Exit(code=1) from None

    # Calculate file size
    file_size_mb = output.stat().st_size / (1024 * 1024)

    out.print(f"  [green]Report saved to {output}[/green]")
    out.print(f"  [dim]File size: {file_size_mb:.1f} MB[/dim]")

    # Summary
    out.print("\n[bold green]Multi-sample report complete![/bold green]")
    out.print("\n[bold]Report contents:[/bold]")
    out.print(f"  - Samples compared: {len(sample_data)}")
    out.print(f"  - Total reads: {total_reads:,}")
    out.print(f"  - Theme: {theme}")
    out.print(f"\n[dim]Open in browser: file://{output.absolute()}[/dim]\n")


@app.command(name="summarize")
def summarize_novel_diversity(
    classifications: Path = typer.Argument(
        ...,
        help="Classification results file (CSV or Parquet)",
        exists=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output file path (TSV or JSON based on --format)",
    ),
    metadata: Path | None = typer.Option(
        None,
        "--metadata", "-m",
        help="Genome metadata file (TSV) for species/genus/family lookups",
        exists=True,
        dir_okay=False,
    ),
    novelty_band_size: float = typer.Option(
        5.0,
        "--band-size", "-b",
        help="Size of novelty bands in percent (default 5.0)",
    ),
    min_cluster_size: int = typer.Option(
        3,
        "--min-cluster-size", "-s",
        help="Minimum reads to form a cluster (default 3)",
    ),
    output_format: str = typer.Option(
        "tsv",
        "--format", "-f",
        help="Output format: 'tsv' or 'json'",
    ),
    json_output: Path | None = typer.Option(
        None,
        "--json", "-j",
        help="Additional JSON output path with full metadata",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Summarize novel diversity into putative taxon clusters.

    Groups reads classified as Novel Species or Novel Genus into clusters
    based on their nearest reference genome and novelty level. This enables
    identification of distinct novel organisms for further investigation.

    Example:

        # Basic summary (TSV output)
        metadarkmatter report summarize classifications.csv \\
            --output novel_summary.tsv \\
            --metadata genome_metadata.tsv

        # With JSON output for programmatic access
        metadarkmatter report summarize classifications.csv \\
            --output novel_summary.tsv \\
            --json novel_summary.json \\
            --metadata genome_metadata.tsv

        # Custom clustering parameters
        metadarkmatter report summarize classifications.csv \\
            --output novel_summary.tsv \\
            --band-size 10.0 \\
            --min-cluster-size 5
    """
    import json

    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Novel Diversity Summary[/bold blue]\n")

    # Validate format
    if output_format not in ("tsv", "json"):
        console.print(
            f"[red]Error: Invalid format '{output_format}'. Use 'tsv' or 'json'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Step 1: Load classification data
    out.print("[bold]Step 1:[/bold] Loading classification data...")

    with spinner_progress("Reading classifications...", console, quiet):
        try:
            df = read_dataframe(classifications)
        except Exception as e:
            console.print(f"[red]Error reading classifications: {e}[/red]")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Loaded {len(df):,} classifications[/green]")

    # Count novel reads
    novel_count = len(df.filter(
        pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
    ))
    out.print(f"  [dim]Novel reads: {novel_count:,}[/dim]")

    if novel_count == 0:
        console.print("[yellow]No novel reads found in classifications.[/yellow]")
        raise typer.Exit(code=0) from None

    # Step 2: Load metadata if provided
    genome_metadata: GenomeMetadata | None = None
    if metadata:
        out.print("\n[bold]Step 2:[/bold] Loading genome metadata...")
        try:
            genome_metadata = GenomeMetadata.from_file(metadata)
            out.print(
                f"  [green]Loaded metadata for {genome_metadata.genome_count} genomes "
                f"({genome_metadata.species_count} species)[/green]"
            )
        except Exception as e:
            console.print(f"[yellow]Warning: Could not load metadata: {e}[/yellow]")

    # Step 3: Cluster novel reads
    out.print("\n[bold]Step 3:[/bold] Clustering novel reads...")

    with spinner_progress("Analyzing novel diversity...", console, quiet):
        try:
            from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer

            analyzer = NovelDiversityAnalyzer(
                classifications=df,
                metadata=genome_metadata,
                novelty_band_size=novelty_band_size,
                min_cluster_size=min_cluster_size,
            )
            clusters = analyzer.cluster_novel_reads()
            summary = analyzer.get_summary()

        except Exception as e:
            console.print(f"[red]Error analyzing novel diversity: {e}[/red]")
            if verbose:
                console.print_exception()
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Found {len(clusters)} clusters[/green]")

    if len(clusters) == 0:
        console.print(
            f"[yellow]No clusters found with >= {min_cluster_size} reads.[/yellow]"
        )
        raise typer.Exit(code=0) from None

    # Show summary
    if verbose:
        out.print(f"\n[bold]Cluster Summary:[/bold]")
        out.print(f"  - Novel Species clusters: {summary.novel_species_clusters}")
        out.print(f"  - Novel Genus clusters: {summary.novel_genus_clusters}")
        out.print(f"  - High confidence: {summary.high_confidence_clusters}")
        out.print(f"  - Medium confidence: {summary.medium_confidence_clusters}")
        out.print(f"  - Low confidence: {summary.low_confidence_clusters}")

    # Step 4: Write output
    out.print("\n[bold]Step 4:[/bold] Writing output...")

    # Ensure output directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    try:
        if output_format == "tsv":
            cluster_df = analyzer.to_dataframe()
            cluster_df.write_csv(output, separator="\t")
            out.print(f"  [green]Saved TSV to {output}[/green]")
        else:
            output_data = analyzer.to_dict()
            output.write_text(json.dumps(output_data, indent=2))
            out.print(f"  [green]Saved JSON to {output}[/green]")

        # Write additional JSON if requested
        if json_output:
            json_output.parent.mkdir(parents=True, exist_ok=True)
            output_data = analyzer.to_dict()
            json_output.write_text(json.dumps(output_data, indent=2))
            out.print(f"  [green]Saved full JSON to {json_output}[/green]")

    except Exception as e:
        console.print(f"[red]Error writing output: {e}[/red]")
        raise typer.Exit(code=1) from None

    # Summary
    out.print("\n[bold green]Novel diversity summary complete![/bold green]")
    out.print("\n[bold]Summary:[/bold]")
    out.print(f"  - Total clusters: {len(clusters)}")
    out.print(f"  - Novel Species: {summary.novel_species_clusters} clusters ({summary.novel_species_reads:,} reads)")
    out.print(f"  - Novel Genus: {summary.novel_genus_clusters} clusters ({summary.novel_genus_reads:,} reads)")
    out.print(f"  - High confidence: {summary.high_confidence_clusters}")
    out.print(f"  - Genera with novel species: {summary.genera_with_novel_species}")
    out.print(f"  - Families with novel genera: {summary.families_with_novel_genera}")
    out.print(f"\n[dim]Output: {output}[/dim]\n")
