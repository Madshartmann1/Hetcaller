import sys

def count_sites_with_het(file_path, window_size, min_sites, min_proportion):
    site_count = 0
    het_count = 0
    current_window_start = 0
    current_window_end = window_size - 1
    current_scaffold = None
    consecutive_low_het_count = 0
    consecutive_low_het_dict = {}

    # Read the input file
    with open(file_path, 'r') as file:
        for line in file:
            try:
                scaffold, site, *other_fields = line.strip().split('\t')
                site = int(site)
            except ValueError:
                # Skip lines with invalid site values
                continue

            # Check if the scaffold has changed
            if scaffold != current_scaffold:
                # Reset the counts for the new scaffold
                site_count = 0
                het_count = 0
                current_window_start = 0
                current_window_end = window_size - 1
                current_scaffold = scaffold
                consecutive_low_het_count = 0

            # Check if the site is within the current window
            if site > current_window_end:
                # Calculate the proportion of "HET" sites in the window (avoid division by zero)
                proportion_het = het_count / site_count if site_count > 0 else 0

                # Check if the window has at least min_sites
                if site_count >= min_sites:
                    # Check if the proportion is less than min_proportion
                    if proportion_het < min_proportion:
                        consecutive_low_het_count += 1
                    else:
                        # Add the count to the dictionary based on the number of consecutive low "HET" windows
                        if consecutive_low_het_count > 0:
                            consecutive_low_het_dict[consecutive_low_het_count] = consecutive_low_het_dict.get(consecutive_low_het_count, 0) + 1

                        consecutive_low_het_count = 0

                    # Print the count for the completed window
                    print(f"{current_scaffold} {current_window_start}-{current_window_end}: {site_count} sites, {het_count} sites with 'HET', Proportion: {proportion_het:.7f}")

                # Reset the counts for the new window
                site_count = 0
                het_count = 0
                current_window_start += window_size
                current_window_end += window_size

            # Increment the site count for the window
            site_count += 1

            # Check if the line contains "HET"
            if "HET_trv" in line:
                het_count += 1

    # Calculate the proportion of "HET" sites in the last window (avoid division by zero)
    proportion_het = het_count / site_count if site_count > 0 else 0

    # Print the count for the last window
    print(f"{current_scaffold} {current_window_start}-{current_window_end}: {site_count} sites, {het_count} sites with 'HET', Proportion: {proportion_het:.7f}")

    # Add the count to the dictionary based on the number of consecutive low "HET" windows
    if consecutive_low_het_count > 0:
        consecutive_low_het_dict[consecutive_low_het_count] = consecutive_low_het_dict.get(consecutive_low_het_count, 0)

    # Print the counts for consecutive low "HET" windows
    print("\nConsecutive low 'HET' windows:")
    for num_consecutive, count in consecutive_low_het_dict.items():
        print(f"{num_consecutive} consecutive windows: {count}")


# Check if the input file path and arguments are provided correctly
if len(sys.argv) < 5:
    print("Please provide the input file path, window size, minimum number of sites, and minimum proportion as arguments.")
    sys.exit(1)

# Get the input file path and arguments from command-line arguments
file_path = sys.argv[1]
window_size = int(sys.argv[2])
min_sites = int(sys.argv[3])
min_proportion = float(sys.argv[4])

# Call the function with the provided arguments
count_sites_with_het(file_path, window_size, min_sites, min_proportion)

