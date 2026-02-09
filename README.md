# Hetcaller

Hetcaller is a Snakemake-based pipeline for analyzing heterozygosity and identifying runs of homozygosity (ROH) from BAM files using ANGSD.

## Features

- **Heterozygosity calling** at customizable thresholds
- **ROH detection** with configurable minimum thresholds
- **Transversion-only analysis** to reduce bias from damage patterns
- **Minor allele frequency analysis** with visualizations
- **Batch processing** of multiple BAM files
- **Flexible region specification** (whole genome, chromosomes, or custom regions)

## Requirements

### Software Dependencies
- Python 3.x with packages:
  - `pyyaml`
  - `snakemake`
- [ANGSD](http://www.popgen.dk/angsd/) (tested with v0.935)
- R (tested with R 4.1.2) with packages:
  - `ggplot2`
  - `dplyr`
  - `readr`
  - `patchwork`
  - `RColorBrewer`
- Standard Unix tools: `bash`, `awk`, `grep`, `paste`, `zcat`

### Installation

1. Clone the repository:
```bash
git clone https://github.com/Madshartmann1/Hetcaller.git
cd Hetcaller
```

2. Make hetcall.py executable:
```bash
chmod +x hetcall.py
```

## Quick Start

### Single BAM Analysis
```bash
./hetcall.py -b sample.bam -o sample01
```

### Multiple BAM Files
```bash
# Create a BAM list file (one BAM per line)
./hetcall.py -l bam_list.txt
```

### With Custom Thresholds
```bash
./hetcall.py -b sample.bam -o sample01 -t 0.05 -t 0.10 -t 0.15
```

### Region-specific Analysis
```bash
./hetcall.py -b sample.bam -o sample01 -r "chr1:1-50000000"
```

## Usage

```
hetcall [OPTIONS]

Required Arguments:
  -b, --bam FILE          BAM file for analysis (or use --bam-list)
  -l, --bam-list FILE     File with BAM paths, one per line
  -o, --out-prefix STR    Output prefix (required for single BAM)

Optional Arguments:
  -d, --outdir DIR        Output directory [default: results]
  -s, --scripts DIR       Directory with required scripts [default: scripts]
  -m, --min-depth INT     Minimum depth for ANGSD [default: 10]
  -t, --threshold FLOAT   Het calling threshold (can be repeated) [default: 0.05]
  -T, --threshold-list FILE  File with thresholds, one per line
  -R, --roh-min FLOAT     ROH minimum threshold [default: 0.2]
  -c, --cores INT         Number of cores [default: 8]
  -r, --regions STR       Region string for ANGSD (e.g., 'chr1:1-200000000')
  -f, --rf FILE           Regions file (one scaffold per line)
  -n, --dry-run           Show what would be done
  -h, --help              Show detailed help
```

## Output Files

For each sample, Hetcaller generates:

- `{prefix}.pos.gz` - Genomic positions analyzed
- `{prefix}.counts.gz` - ANGSD allele counts
- `{prefix}.basecall_{threshold}_mincov{depth}.txt` - Base calls at each position
- `{prefix}_het{threshold}_thres{roh_min}.ROH.txt` - ROH analysis (all SNPs)
- `{prefix}_het{threshold}_thres{roh_min}.ROH_trv.txt` - ROH analysis (transversions only)
- `{prefix}.minorfreq*.txt` - Minor allele frequency distributions
- `plots/{prefix}_minorfreq_plots.pdf` - Visualization plots

## BAM List Format

Each line in the BAM list file should contain:
```
/path/to/sample.bam [optional_prefix]
```

If no prefix is provided, the BAM filename (without .bam extension) will be used.

## Examples

See the `examples/` directory for:
- `example_config.yaml` - Sample configuration file
- `example_bam_list.txt` - Sample BAM list format
- `example_thresholds.txt` - Sample threshold file


## Contact

madhar@dtu.dk