# Changelog

All notable changes to metadarkmatter are documented in this file.

## [Unreleased] - 2026-01-20

### Added

**BLAST FASTQ Support**
- BLAST now accepts FASTQ files directly (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`)
- Automatic conversion to FASTA using `seqtk` (primary) or Python fallback
- Transparent temporary file management with automatic cleanup
- Users no longer need to manually convert FASTQ to FASTA

**MMseqs2 Integration**
- MMseqs2 support as BLAST alternative for very large datasets (>100K reads)
- Multi-step workflow with query database caching for performance
- BLAST-compatible output format for seamless pipeline integration
- CLI commands: `metadarkmatter mmseqs2 makedb` and `metadarkmatter mmseqs2 search`

**Project Organization**
- External software dependencies clearly documented in all main files
- Organized test data into `test_data/` directory
- Added `.gitignore` for Python artifacts and test outputs

### Changed

**Documentation Updates**
- Streamlined MMseqs2 documentation with realistic performance expectations
- MMseqs2 is 20x slower than BLAST for small datasets (<10K reads)
- MMseqs2 provides 5-15x speedup only for datasets >100K reads
- Updated all documentation to reflect BLAST as primary method
- Made `seqtk` optional (only needed for assembly workflows)
- Reorganized documentation with tutorial as primary entry point

**Performance Guidance**
- Added realistic performance tables based on actual testing
- Clear dataset size thresholds for tool selection
- BLAST recommended for <100K reads, MMseqs2 for larger datasets

### Fixed

**MMseqs2 Issues**
- Fixed duplicate `-q` parameter warning in CLI
- Implemented query database caching to prevent 30+ min overhead per run
- Added missing `--search-type` parameter to search command
- Fixed `get_executable_path()` method name error

**Documentation Accuracy**
- Corrected misleading "100-1000x speedup" claims for all datasets
- Removed overstated MMseqs2 performance claims
- Added warnings about dataset size requirements

### Developer Notes

**File Reorganization**
- Moved documentation summaries from root to `docs/`
- Consolidated test files into `test_data/demo_outputs/` and `test_data/e2e_tests/`
- Created CHANGELOG.md (this file) from CLEANUP_SUMMARY.md, DOCUMENTATION_UPDATE_SUMMARY.md

**External Dependencies**
- Required: kraken2, krakentools, blast, skani
- Optional: mmseqs2 (>100K reads), seqtk (assembly workflows)

## Implementation Details

For detailed implementation notes, see:
- `docs/MMSEQS2_IMPLEMENTATION_SUMMARY.md` - Complete MMseqs2 integration details
- `docs/DOCUMENTATION_UPDATE_SUMMARY.md` - Documentation update details
