#!/usr/bin/env python3
# hetcall: CLI wrapper that runs Snakemake with your config

import argparse, subprocess, sys, yaml, os, tempfile, shutil, shlex

def show_help():
    """Show detailed help with examples"""
    help_text = """
hetcall - Run ANGSD→basecalls→ROH analysis via Snakemake

USAGE:
    hetcall [OPTIONS]

REQUIRED ARGUMENTS:
    -b, --bam FILE          BAM file for analysis (or use --bam)
    -l, --bam-list FILE     File with BAM paths, one per line (or use --bam-list)
    -o, --out-prefix STR    Output prefix for all files (ignored if using --bam-list)

OPTIONAL ARGUMENTS:
    -d, --outdir DIR        Output directory [default: results]
    -s, --scripts DIR       Directory with required scripts [default: diversity]
    -m, --min-depth INT     Minimum depth for ANGSD [default: 10]
    -t, --threshold FLOAT   Het calling threshold (can be repeated) [default: 0.05]
    -T, --threshold-list FILE  File with thresholds, one per line
    -R, --roh-min FLOAT     ROH minimum threshold (can be repeated) [default: 0.2]
    -L, --roh-list FILE     File with ROH thresholds, one per line
    -c, --cores INT         Number of cores for Snakemake [default: 8]
    
REGION SPECIFICATION (choose one):
    -r, --regions STR       Region string for ANGSD -r (e.g., 'chr1:' or 'chr1:1-200000000')
    -f, --rf FILE           Regions file for ANGSD -rf (one scaffold per line)

WORKFLOW OPTIONS:
    -S, --snakefile FILE    Snakefile path [default: Snakefile]
    -C, --configfile FILE  Config file path [default: config.yaml in current dir]
    -n, --dry-run          Show what would be done without executing
    -u, --unlock           Unlock working directory
    -F, --force            Overwrite existing config file
    -h, --help             Show this help message

EXAMPLES:
    # Basic usage with single BAM
    hetcall -b sample.bam -o sample01
    
    # Multiple BAMs from list (runs sequentially)
    hetcall -l bam_list.txt
    
    # Multiple thresholds
    hetcall -b sample.bam -o sample01 -t 0.05 -t 0.1 -t 0.15
    
    # Thresholds from file
    hetcall -b sample.bam -o sample01 -T thresholds.txt
    
    # Multiple ROH thresholds
    hetcall -b sample.bam -o sample01 -R 0.1 -R 0.2 -R 0.3
    
    # Use specific chromosome region
    hetcall -b sample.bam -o sample01 -r "chr1:1-50000000"
    
    # Dry run with multiple BAMs
    hetcall -l bam_list.txt -n

BAM LIST FORMAT:
    Each line: /path/to/file.bam [optional_prefix]
    If no prefix provided, uses basename without .bam extension
    
    Example bam_list.txt:
    /data/sample1.bam sample1
    /data/sample2.bam sample2  
    /data/sample3.bam
"""
    print(help_text)

def read_file_list(filepath, list_type="values"):
    """Read a file with one item per line, return list"""
    if not os.path.exists(filepath):
        print(f"Error: {list_type} file {filepath} not found", file=sys.stderr)
        sys.exit(1)
    
    items = []
    with open(filepath) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line and not line.startswith('#'):
                if list_type == "bams":
                    # Handle BAM list format: /path/to/file.bam [optional_prefix]
                    parts = line.split()
                    bam_path = parts[0]
                    prefix = parts[1] if len(parts) > 1 else os.path.splitext(os.path.basename(bam_path))[0]
                    items.append((bam_path, prefix))
                else:
                    # Simple list of values
                    try:
                        items.append(float(line))
                    except ValueError:
                        print(f"Error: Invalid {list_type} value '{line}' on line {line_num} in {filepath}", file=sys.stderr)
                        sys.exit(1)
    
    if not items:
        print(f"Error: No valid entries found in {filepath}", file=sys.stderr)
        sys.exit(1)
    
    return items

def run_single_bam(bam_path, prefix, args, base_cfg):
    """Run analysis for a single BAM file"""
    print(f"\n{'='*60}")
    print(f"Processing: {bam_path} -> {prefix}")
    print(f"{'='*60}")
    
    # Create config for this BAM
    cfg = base_cfg.copy()
    cfg["bam"] = bam_path
    cfg["prefix"] = prefix
    
    # Config file path
    cfg_path = args.configfile or f"config_{prefix}.yaml"
    
    # Check if config exists and handle --force
    if os.path.exists(cfg_path) and not args.force:
        print(f"Warning: Config file {cfg_path} already exists. Use --force to overwrite.", file=sys.stderr)
        response = input("Continue anyway? [y/N]: ").lower().strip()
        if response not in ['y', 'yes']:
            print("Skipping this BAM.")
            return False
    
    # Write config
    with open(cfg_path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    
    print(f"Config written to: {cfg_path}")
    
    # Build snakemake command
    cmd = [
        "snakemake",
        "-s", args.snakefile,
        "--cores", str(args.cores),
        "--rerun-incomplete",
        "--printshellcmds",
        "--configfile", cfg_path,
    ]
    
    if args.dry_run:
        cmd.append("--dry-run")
    
    if args.unlock:
        cmd.append("--unlock")
    
    print(f"Running: {' '.join(cmd)}")
    
    if args.dry_run:
        print("(Dry run - no actual execution)")
        return True
    
    # Run the command
    result = subprocess.call(cmd)
    if result != 0:
        print(f"Error: Snakemake failed for {prefix} (exit code {result})", file=sys.stderr)
        return False
    
    print(f"Successfully completed: {prefix}")
    return True

def main():
    # Custom help handling
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ['-h', '--help']):
        show_help()
        sys.exit(0)
    
    p = argparse.ArgumentParser(
        prog="hetcall", 
        description="Run ANGSD→basecalls→ROH analysis via Snakemake",
        add_help=False  # We handle help ourselves
    )
    
    # Required arguments (now mutually exclusive)
    required = p.add_argument_group('Required arguments (choose one)')
    bam_group = required.add_mutually_exclusive_group(required=True)
    bam_group.add_argument("-b", "--bam", metavar="FILE",
                          help="BAM file for analysis")
    bam_group.add_argument("-l", "--bam-list", metavar="FILE",
                          help="File with BAM paths, one per line")
    
    required.add_argument("-o", "--out-prefix", metavar="STR",
                         help="Output prefix (required for single BAM, ignored for BAM list)")
    
    # Optional arguments
    optional = p.add_argument_group('Optional arguments')
    optional.add_argument("-d", "--outdir", default="results", metavar="DIR",
                         help="Output directory (default: %(default)s)")
    optional.add_argument("-s", "--scripts", default="diversity", metavar="DIR",
                         help="Directory with required scripts (default: %(default)s)")
    optional.add_argument("-m", "--min-depth", type=int, default=10, metavar="INT",
                         help="Minimum depth for ANGSD (default: %(default)s)")
    optional.add_argument("-c", "--cores", type=int, default=8, metavar="INT",
                         help="Number of cores for Snakemake (default: %(default)s)")
    
    # Threshold options
    thresh_group = optional.add_mutually_exclusive_group()
    thresh_group.add_argument("-t", "--threshold", type=float, action="append", 
                             default=[], metavar="FLOAT",
                             help="Het calling threshold (can be repeated, default: 0.05)")
    thresh_group.add_argument("-T", "--threshold-list", metavar="FILE",
                             help="File with thresholds, one per line")
    
    # ROH threshold options  
    roh_group = optional.add_mutually_exclusive_group()
    roh_group.add_argument("-R", "--roh-min", type=float, action="append",
                          default=[], metavar="FLOAT",
                          help="ROH minimum threshold (can be repeated, default: 0.2)")
    roh_group.add_argument("-L", "--roh-list", metavar="FILE",
                          help="File with ROH thresholds, one per line")
    
    # Region specification
    regions = p.add_argument_group('Region specification (choose one)')
    region_group = regions.add_mutually_exclusive_group()
    region_group.add_argument("-r", "--regions", metavar="STR",
                             help="Region string for ANGSD -r (e.g., 'chr1:' or 'chr1:1-200000000')")
    region_group.add_argument("-f", "--rf", metavar="FILE",
                             help="Regions file for ANGSD -rf (one scaffold per line)")
    
    # Workflow options
    workflow = p.add_argument_group('Workflow options')
    workflow.add_argument("-S", "--snakefile", default="Snakefile", metavar="FILE",
                         help="Snakefile path (default: %(default)s)")
    workflow.add_argument("-C", "--configfile", metavar="FILE",
                         help="Config file path (default: auto-generated)")
    workflow.add_argument("-n", "--dry-run", action="store_true",
                         help="Show what would be done without executing")
    workflow.add_argument("-u", "--unlock", action="store_true",
                         help="Unlock working directory")
    workflow.add_argument("-F", "--force", action="store_true",
                         help="Overwrite existing config file")
    workflow.add_argument("-h", "--help", action="store_true",
                         help="Show detailed help message")
    
    args = p.parse_args()
    
    # Handle help
    if args.help:
        show_help()
        sys.exit(0)
    
    # Validate single BAM requires prefix
    if args.bam and not args.out_prefix:
        print("Error: --out-prefix is required when using --bam", file=sys.stderr)
        sys.exit(1)
    
    # Handle thresholds
    if args.threshold_list:
        thresholds = read_file_list(args.threshold_list, "thresholds")
    elif args.threshold:
        thresholds = args.threshold
    else:
        thresholds = [0.05]
    
    # Handle ROH thresholds
    if args.roh_list:
        roh_mins = read_file_list(args.roh_list, "ROH thresholds")
    elif args.roh_min:
        roh_mins = args.roh_min
    else:
        roh_mins = [0.2]
    
    # Base config template
    base_cfg = {
        "rf": args.rf,
        "regions": args.regions,
        "scripts": args.scripts,
        "outdir": args.outdir,
        "min_depth": args.min_depth,
        "thresholds": thresholds,
        "roh_min": roh_mins[0] if len(roh_mins) == 1 else roh_mins,  # Single value or list
        "angsd_threads": max(1, min(args.cores, 64)),
    }
    
    # Handle BAMs
    if args.bam_list:
        bam_list = read_file_list(args.bam_list, "bams")
        print(f"Found {len(bam_list)} BAM files to process")
        
        success_count = 0
        for bam_path, prefix in bam_list:
            if run_single_bam(bam_path, prefix, args, base_cfg):
                success_count += 1
        
        print(f"\n{'='*60}")
        print(f"Completed {success_count}/{len(bam_list)} BAM files successfully")
        if success_count < len(bam_list):
            sys.exit(1)
            
    else:
        # Single BAM
        success = run_single_bam(args.bam, args.out_prefix, args, base_cfg)
        if not success:
            sys.exit(1)

if __name__ == "__main__":
    main()
