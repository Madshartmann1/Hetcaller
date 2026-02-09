# GitHub Setup Instructions for Hetcaller

Your Hetcaller repository is ready! Here's how to push it to GitHub:

## Step 1: Create GitHub Repository

1. Go to https://github.com/new
2. Repository name: `Hetcaller`
3. Description: "Snakemake pipeline for heterozygosity and ROH analysis from BAM files using ANGSD"
4. Keep it Public or Private (your choice)
5. **DO NOT** initialize with README, .gitignore, or license (we already have these)
6. Click "Create repository"

## Step 2: Push to GitHub

After creating the repository, run these commands:

```bash
cd /home/projects2/Biodiversity_Extinction/Hetcaller

# Add the remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/Hetcaller.git

# Rename branch to main (optional, GitHub default)
git branch -M main

# Push to GitHub
git push -u origin main
```

## Alternative: Using SSH

If you have SSH keys set up:

```bash
cd /home/projects2/Biodiversity_Extinction/Hetcaller
git remote add origin git@github.com:YOUR_USERNAME/Hetcaller.git
git branch -M main
git push -u origin main
```

## What's Included

The repository contains:

```
Hetcaller/
├── .gitignore               # Excludes results, data files, etc.
├── README.md                # Main documentation
├── hetcall.py              # Main CLI wrapper script
├── Snakefile               # Snakemake workflow definition
├── examples/               
│   ├── example_config.yaml      # Example configuration
│   ├── example_bam_list.txt     # Example BAM list format
│   └── example_thresholds.txt   # Example threshold file
└── scripts/                
    ├── Basecall.sh         # Base calling script
    ├── header.txt          # Header for output files
    ├── ROH.py              # ROH detection (all SNPs)
    ├── ROH_trv.py          # ROH detection (transversions only)
    └── plot_minorfreq.R    # Visualization script
```

## Next Steps (Optional)

1. **Add a LICENSE**: Consider adding MIT, GPL, or another open-source license
2. **Add topics**: On GitHub, add topics like `genomics`, `heterozygosity`, `roh`, `snakemake`, `angsd`
3. **Create releases**: Tag versions as you make updates
4. **Enable Issues**: Allow users to report bugs and request features

## Verification

After pushing, verify your repository at:
https://github.com/YOUR_USERNAME/Hetcaller

The README.md will be displayed on the main page with installation and usage instructions.
