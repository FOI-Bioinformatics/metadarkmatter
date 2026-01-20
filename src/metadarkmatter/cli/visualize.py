"""
Visualize command for creating recruitment plots.

Creates interactive Plotly visualizations from BAM alignment data
showing read recruitment patterns across reference genomes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import typer
from rich.console import Console

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.core.io_utils import read_dataframe
from metadarkmatter.core.recruitment import (
    aggregate_by_genome,
    load_recruitment_data,
)
from metadarkmatter.external import Samtools

app = typer.Typer(
    name="visualize",
    help="Create recruitment plots from BAM alignment data",
    no_args_is_help=True,
)

console = Console()


OutputFormat = Literal["html", "png", "json"]


@app.command(name="recruitment")
def recruitment_plot(
    bam: Path = typer.Option(
        ...,
        "--bam", "-b",
        help="Input BAM file from competitive mapping",
        exists=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output plot file path",
    ),
    output_format: str = typer.Option(
        "html",
        "--format", "-f",
        help="Output format: 'html' (interactive), 'png' (static), or 'json' (data)",
    ),
    genome: list[str] | None = typer.Option(
        None,
        "--genome", "-g",
        help="Specific genome(s) to plot (can be repeated)",
    ),
    min_identity: float = typer.Option(
        70.0,
        "--min-identity",
        help="Minimum percent identity to include",
        min=0.0,
        max=100.0,
    ),
    min_mapq: int = typer.Option(
        0,
        "--min-mapq",
        help="Minimum mapping quality to include",
        min=0,
    ),
    show_bands: bool = typer.Option(
        True,
        "--show-bands/--hide-bands",
        help="Show identity threshold bands",
    ),
    known_threshold: float = typer.Option(
        98.0,
        "--known-threshold",
        help="Identity threshold for known species",
        min=0.0,
        max=100.0,
    ),
    novel_species_threshold: float = typer.Option(
        85.0,
        "--novel-species-threshold",
        help="Identity threshold for novel species",
        min=0.0,
        max=100.0,
    ),
    novel_genus_threshold: float = typer.Option(
        75.0,
        "--novel-genus-threshold",
        help="Identity threshold for novel genus",
        min=0.0,
        max=100.0,
    ),
    max_points: int = typer.Option(
        100000,
        "--max-points",
        help="Maximum points to plot (subsamples if exceeded)",
        min=1000,
    ),
    width: int = typer.Option(
        1200,
        "--width",
        help="Plot width in pixels",
        min=400,
    ),
    height: int = typer.Option(
        800,
        "--height",
        help="Plot height in pixels",
        min=300,
    ),
    export_tsv: Path | None = typer.Option(
        None,
        "--export-tsv",
        help="Export data as TSV for anvi'o compatibility",
    ),
    multi_panel: bool = typer.Option(
        False,
        "--multi-panel",
        help="Create multi-panel plot with one subplot per genome",
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
    Create recruitment plot from BAM alignment data.

    This command generates an interactive visualization showing read
    recruitment patterns across reference genomes. The plot displays
    percent identity vs genome position, with horizontal bands
    indicating taxonomic classification thresholds.

    The visualization helps identify:
    - Known species (>98% identity band)
    - Novel species (85-98% identity band)
    - Novel genera (75-85% identity band)

    Example:

        # Create interactive HTML plot
        metadarkmatter visualize recruitment \\
            --bam sample_mapped.bam \\
            --output recruitment.html

        # Create PNG for a specific genome
        metadarkmatter visualize recruitment \\
            --bam sample_mapped.bam \\
            --output recruitment.png \\
            --format png \\
            --genome GCF_000123456.1

        # Export data for anvi'o
        metadarkmatter visualize recruitment \\
            --bam sample_mapped.bam \\
            --output recruitment.html \\
            --export-tsv recruitment_data.tsv
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Recruitment Plot[/bold blue]\n")

    # Validate format
    output_format_lower = output_format.lower()
    if output_format_lower not in ("html", "png", "json"):
        console.print(
            f"[red]Error: Invalid format '{output_format}'. "
            f"Use 'html', 'png', or 'json'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Check samtools availability
    if not Samtools.check_available():
        console.print(
            "[red]Error: samtools not found in PATH[/red]\n"
            "[dim]Install with: conda install -c bioconda samtools[/dim]"
        )
        raise typer.Exit(code=1) from None

    # Try to import plotly
    try:
        from metadarkmatter.visualization import (
            IdentityBand,
            RecruitmentPlotGenerator,
        )
    except ImportError:
        console.print(
            "[red]Error: Plotly is required for visualization[/red]\n"
            "[dim]Install with: pip install plotly[/dim]"
        )
        raise typer.Exit(code=1) from None

    # Ensure output directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    # Load recruitment data from BAM
    out.print("[bold]Step 1:[/bold] Loading alignment data...")

    with spinner_progress("Reading BAM file...", console, quiet):
        try:
            data = load_recruitment_data(
                bam_path=bam,
                min_mapq=min_mapq,
                min_identity=min_identity,
                max_records=max_points * 10,  # Load more than needed for filtering
            )
        except Exception as e:
            console.print(f"\n[red]Error reading BAM file: {e}[/red]")
            raise typer.Exit(code=1) from None

    num_alignments = len(data)
    out.print(f"  [green]Loaded {num_alignments:,} alignments[/green]")

    if num_alignments == 0:
        console.print("[yellow]Warning: No alignments found matching filters[/yellow]")
        raise typer.Exit(code=0)

    # Show genome statistics
    if verbose:
        genome_stats = aggregate_by_genome(data)
        out.print("\n[dim]Top genomes by alignment count:[/dim]")
        for row in genome_stats.head(5).iter_rows(named=True):
            out.print(
                f"  [dim]{row['genome_name']}: {row['num_reads']:,} reads, "
                f"mean identity {row['mean_identity']:.1f}%[/dim]"
            )

    # Filter to specific genomes if requested
    if genome:
        data = data.filter(data["genome_name"].is_in(genome))
        out.print(f"  [dim]Filtered to {len(data):,} alignments for {len(genome)} genome(s)[/dim]")

    # Create custom bands if thresholds differ from defaults
    custom_bands = (
        IdentityBand(
            name="Known Species",
            min_identity=known_threshold,
            max_identity=100.0,
            color="#2ecc71",
        ),
        IdentityBand(
            name="Novel Species",
            min_identity=novel_species_threshold,
            max_identity=known_threshold,
            color="#f39c12",
        ),
        IdentityBand(
            name="Novel Genus",
            min_identity=novel_genus_threshold,
            max_identity=novel_species_threshold,
            color="#e74c3c",
        ),
    )

    # Create plot
    out.print("\n[bold]Step 2:[/bold] Generating plot...")

    with spinner_progress("Creating visualization...", console, quiet):
        generator = RecruitmentPlotGenerator(data=data, bands=custom_bands)

        try:
            if multi_panel:
                fig = generator.create_multi_genome_figure(
                    genomes=list(genome) if genome else None,
                    show_bands=show_bands,
                    max_points_per_genome=max_points // 5,
                )
            else:
                fig = generator.create_figure(
                    show_bands=show_bands,
                    max_points=max_points,
                    width=width,
                    height=height,
                )

            # Save plot
            if output_format_lower == "html":
                fig.write_html(output, include_plotlyjs=True)
            elif output_format_lower == "png":
                fig.write_image(output)
            elif output_format_lower == "json":
                fig.write_json(output)

        except Exception as e:
            console.print(f"\n[red]Error creating plot: {e}[/red]")
            if verbose:
                console.print_exception()
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Plot saved to {output}[/green]")

    # Export TSV for anvi'o if requested
    if export_tsv:
        out.print("\n[bold]Step 3:[/bold] Exporting data for anvi'o...")

        with spinner_progress("Writing TSV...", console, quiet):
            try:
                export_tsv.parent.mkdir(parents=True, exist_ok=True)
                generator.export_for_anvio(export_tsv)
            except Exception as e:
                console.print(f"\n[yellow]Warning: Failed to export TSV: {e}[/yellow]")

        out.print(f"  [green]TSV exported to {export_tsv}[/green]")

    # Summary
    out.print("\n[bold green]Visualization complete![/bold green]")
    out.print("\n[bold]Output files:[/bold]")
    out.print(f"  Plot: {output}")
    if export_tsv:
        out.print(f"  Data: {export_tsv}")
    out.print()


@app.command(name="summary")
def summary_plot(
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
        help="Output plot file path",
    ),
    output_format: str = typer.Option(
        "html",
        "--format", "-f",
        help="Output format: 'html' or 'png'",
    ),
    width: int = typer.Option(
        1000,
        "--width",
        help="Plot width in pixels",
    ),
    height: int = typer.Option(
        600,
        "--height",
        help="Plot height in pixels",
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
    Create summary visualization from classification results.

    Generates bar charts and pie charts showing the distribution
    of taxonomic classifications.

    Example:

        metadarkmatter visualize summary \\
            --classifications results.csv \\
            --output summary.html
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Classification Summary[/bold blue]\n")

    # Try to import plotly
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        console.print(
            "[red]Error: Plotly is required for visualization[/red]\n"
            "[dim]Install with: pip install plotly[/dim]"
        )
        raise typer.Exit(code=1) from None

    # Load classification data
    out.print("[bold]Step 1:[/bold] Loading classification data...")

    try:
        df = read_dataframe(classifications)
    except Exception as e:
        console.print(f"[red]Error reading file: {e}[/red]")
        raise typer.Exit(code=1) from None

    out.print(f"  [green]Loaded {len(df):,} classifications[/green]")

    # Calculate summary statistics
    counts = df.group_by("taxonomic_call").len().sort("len", descending=True)
    categories = counts["taxonomic_call"].to_list()
    values = counts["len"].to_list()
    total = sum(values)

    # Color mapping
    colors = {
        "Known Species": "#2ecc71",
        "Novel Species": "#f39c12",
        "Novel Genus": "#e74c3c",
        "Conserved Region": "#95a5a6",
    }

    # Create figure with subplots
    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "bar"}, {"type": "pie"}]],
        subplot_titles=("Classification Counts", "Classification Distribution"),
    )

    # Bar chart
    fig.add_trace(
        go.Bar(
            x=categories,
            y=values,
            marker_color=[colors.get(c, "#333333") for c in categories],
            text=[f"{v:,} ({100*v/total:.1f}%)" for v in values],
            textposition="outside",
        ),
        row=1,
        col=1,
    )

    # Pie chart
    fig.add_trace(
        go.Pie(
            labels=categories,
            values=values,
            marker_colors=[colors.get(c, "#333333") for c in categories],
            textinfo="percent+label",
            hole=0.3,
        ),
        row=1,
        col=2,
    )

    # Update layout
    fig.update_layout(
        title={
            "text": "Taxonomic Classification Summary",
            "x": 0.5,
            "xanchor": "center",
        },
        width=width,
        height=height,
        template="plotly_white",
        showlegend=False,
    )

    # Save
    output.parent.mkdir(parents=True, exist_ok=True)

    output_format_lower = output_format.lower()
    if output_format_lower == "html":
        fig.write_html(output, include_plotlyjs=True)
    elif output_format_lower == "png":
        fig.write_image(output)
    else:
        console.print(f"[red]Error: Invalid format '{output_format}'[/red]")
        raise typer.Exit(code=1) from None

    out.print(f"\n[bold green]Summary plot saved to {output}[/bold green]\n")
