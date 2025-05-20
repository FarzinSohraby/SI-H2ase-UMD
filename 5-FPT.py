import os
import pandas as pd

# Define directories and file structure
base_dir = "/media/farzin/Elements/Hyd-Binding/Hyd-WT-100H2"  # Change this to your actual base directory
replica_dirs = [f"Replica{i}" for i in range(1, 76)]
file_name = "H2_distances.csv"

# Define thresholds
lower_threshold = 5
upper_threshold = 50

# Initialize a report list to store results from all replicas
full_report = []

# Process each replica directory
for replica in replica_dirs:
    # Construct the file path
    file_path = os.path.join(base_dir, replica, file_name)
    
    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue
    
    # Load the CSV file
    data = pd.read_csv(file_path)
    
    # Process each ligand column (ignoring the "Frame" column)
    for ligand in data.columns[1:]:
        distances = data[ligand]
        
        # Find the first frame above 50 Å
        start_frame = (distances > upper_threshold).idxmax()
        
        # Find the first frame where the distance goes below 5 Å
        crossing_frame = (distances <= lower_threshold).idxmax()
        
        # If crossing happens after the start frame, proceed
        if start_frame < crossing_frame:
            # Calculate the number of frames the ligand stays below 5 Å
            below_5_frames = (distances[crossing_frame:] <= lower_threshold).sum()
            
            # Append the results to the report
            full_report.append({
                "Replica": replica,
                "Ligand": ligand,
                "Start Frame": start_frame,
                "Crossing Frame": crossing_frame,
                "Time Below 5 Å (frames)": below_5_frames,
                "Start Distance (Å)": distances.iloc[start_frame],
                "Crossing Distance (Å)": distances.iloc[crossing_frame],
            })

# Convert the full report into a DataFrame for better readability
full_report_df = pd.DataFrame(full_report)

# Sort the results for easier analysis
full_report_df.sort_values(by=["Replica", "Ligand"], inplace=True)
full_report_df.reset_index(drop=True, inplace=True)

# Save the final report to a CSV file for further analysis
output_path = os.path.join(base_dir, "ligand_analysis_report.csv")
full_report_df.to_csv(output_path, index=False)

print(f"Analysis complete. Report saved to {output_path}")

