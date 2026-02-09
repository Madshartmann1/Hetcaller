# Requirements for Hetcaller

## Required Files

The following files are essential for Hetcaller to run:

### Core Scripts
- `hetcall.py` - Main CLI wrapper
- `Snakefile` - Snakemake workflow definition

### Supporting Scripts (in scripts/ directory)
- `Basecall.sh` - Converts allele counts to base calls
- `header.txt` - Column headers for output
- `ROH.py` - ROH detection using all SNPs
- `ROH_trv.py` - ROH detection using transversions only
- `plot_minorfreq.R` - Creates visualization plots

## Software Dependencies

### Python 3.x
Required Python packages:
```bash
pip install pyyaml snakemake
```

Or with conda:
```bash
conda install -c conda-forge -c bioconda snakemake pyyaml
```

### ANGSD (v0.935 or later)
Download and install from: http://www.popgen.dk/angsd/

The Snakefile currently expects ANGSD at:
- `/home/ctools/angsd-0.935/angsd`

You can modify this path in the Snakefile or install ANGSD to this location.

### R (4.1.2 or later)
Required R packages:
```R
install.packages(c("ggplot2", "dplyr", "readr", "patchwork", "RColorBrewer"))
```

The Snakefile looks for Rscript at:
- `/home/ctools/R-4.1.2/bin/Rscript` (preferred)
- `Rscript` (system default)

### System Tools
Standard Unix utilities (typically pre-installed):
- bash
- awk
- grep
- paste
- zcat
- tail
- wc
- sort
- uniq

## Input Requirements

### BAM File
- Aligned sequencing reads in BAM format
- Should be indexed (`.bai` file)
- Quality filtering recommended before analysis

### Reference Genome (optional)
- Required if using region specification (`-r` or `-rf` options)
- Should match the reference used for BAM alignment

## Hardware Requirements

- **Memory**: Depends on genome size and coverage
  - Minimum: 8 GB RAM
  - Recommended: 16+ GB RAM for large genomes
- **CPU**: Multi-threaded (default 8 cores, configurable with `-c`)
- **Disk Space**: 
  - Input BAM size × 2-3 for intermediate files
  - Additional space for results

## Testing Installation

After installing dependencies, test the setup:

```bash
# Check Python dependencies
python3 -c "import yaml, snakemake; print('Python dependencies OK')"

# Check ANGSD
/home/ctools/angsd-0.935/angsd --version

# Check R packages
Rscript -e "library(ggplot2); library(dplyr); library(readr); library(patchwork); library(RColorBrewer); cat('R packages OK\n')"

# Show hetcall help
./hetcall.py --help
```

## Common Issues

1. **ANGSD not found**: Update the path in Snakefile line ~100
2. **Rscript not found**: Install R or update path in Snakefile
3. **Python module errors**: Install missing packages with pip or conda
4. **Permission denied**: Make scripts executable with `chmod +x hetcall.py scripts/*.sh scripts/*.R`
