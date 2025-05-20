import MDAnalysis as mda
import numpy as np
import os
import csv

# Base directory containing all replica directories
base_dir = "/media/farzin/Elements/NewHyd-WT-100H2"

# Summary CSV output file
summary_csv = os.path.join(base_dir, "summary_distances.csv")

# Create or empty the summary CSV file
with open(summary_csv, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(["Replica", "H2 Molecule", "Min Distance (A)", "Frame Number"])

# Loop through each replica directory (Replica1 to Replica75)
for replica_num in range(1, 76):
    replica_dir = os.path.join(base_dir, f"Replica{replica_num}")
    
    # Set file paths specific to each replica
    xtc_file = os.path.join(replica_dir, "NewHyd-WT-H2-Bind-250ns-pbc-ligs.xtc")
    tpr_file = "/media/farzin/Elements/NewHyd-WT-100H2/NewHyd-WT-H2-Bind-250ns-ligs.gro"
    
    # Load the universe (trajectory and topology)
    u = mda.Universe(tpr_file, xtc_file)

    # Select MNI group (assuming it exists in your index.ndx file)
    MNI = u.select_atoms("resname MNI")  # Modify this if needed based on your system

    # Initialize variables to track the minimum distance, corresponding molecule, and frame number
    min_distance = np.inf
    min_molecule = None
    min_frame = None

    # CSV file to write individual H2 distances for the current replica
    replica_csv = os.path.join(replica_dir, "H2_distances.csv")

    # Prepare headers for the H2 molecules (columns)
    H2_molecules = [f"H2_{resid}" for resid in range(799, 899)]
    headers = ["Frame"] + H2_molecules

    # Dictionary to store the distances for each frame
    frame_distances = {}

    # Loop through H2 molecules, which are residues 816 to 915
    for resid in range(799, 899):
        H2 = u.select_atoms(f"resid {resid}")

        # Calculate COM distances for each frame
        for ts in u.trajectory:
            com_H2 = H2.center_of_mass()
            com_MNI = MNI.center_of_mass()
            distance = np.linalg.norm(com_H2 - com_MNI)

            # Store the distance in the dictionary for this frame
            if ts.frame not in frame_distances:
                frame_distances[ts.frame] = [distance]
            else:
                frame_distances[ts.frame].append(distance)

            # Update the minimum distance, corresponding molecule, and frame number
            if distance < min_distance:
                min_distance = distance
                min_molecule = f"H2_{resid}"  # Using residue ID for naming
                min_frame = ts.frame

    # Write distances to the replica CSV in the correct format (each H2 molecule in a separate column)
    with open(replica_csv, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(headers)  # Write header row

        # Write the distances for each frame
        for frame, distances in frame_distances.items():
            csv_writer.writerow([frame] + distances)

    # Append the result for this replica to the summary CSV
    with open(summary_csv, mode='a', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow([f"Replica{replica_num}", min_molecule, min_distance, min_frame])

print(f"Results have been written to {summary_csv}")

