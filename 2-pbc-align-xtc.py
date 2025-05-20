import os
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align

# Paths
csv_path = '/media/farzin/Elements/NewHyd-WT-100H2-R2/ligand_analysis_report-NoT10-6.csv'
replica_dir = '/media/farzin/Elements/NewHyd-WT-100H2-R2'
output_dir = '/media/farzin/Elements/NewHyd-WT-100H2-R2/Aligned-xtc-tauramd-like'

# Load CSV data
csv_data = pd.read_csv(csv_path)

# Set up directories for saving processed trajectories
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Get unique replica numbers from the CSV
unique_replicas = csv_data['Replica'].unique()

# Initialize reference Universe (first replica)
ref_universe = None
ref_atoms = None

# Iterate over each unique replica and process trajectories
for replica_number in unique_replicas:
    replica_name = f'Replica{replica_number}'  # Match exact format in the CSV
    traj_path = os.path.join(replica_dir, replica_name, 'NewHyd-WT-H2-Bind-250ns-pbc-ligs.xtc')  # New trajectory
    top_path = '/media/farzin/Elements/NewHyd-WT-100H2-R2/NewHyd-WT-H2-Bind-250ns-ligs.gro'  # Topology file

    # Check if the H2 molecule data is available for this replica
    h2_data = csv_data[csv_data['Replica'] == replica_number]
    if h2_data.empty:
        print(f"No H2 molecule data found for {replica_name}. Skipping this replica.")
        continue

    # Iterate over each H2 molecule in the replica (if there are multiple)
    for _, row in h2_data.iterrows():
        h2_molecule = int(row['H2 Molecule'])  # Convert to integer

        # Load the trajectory
        try:
            u = mda.Universe(top_path, traj_path)
        except Exception as e:
            print(f"Error loading trajectory or topology for {replica_name}: {e}")
            continue

        # Select protein and the specific H2 molecule by residue ID
        protein = u.select_atoms('resid 1:798')
        h2_molecule_atoms = u.select_atoms(f"resid {h2_molecule}")  # This will now be an integer

        # Check if the selection was successful
        if len(h2_molecule_atoms) == 0:
            print(f"H2 molecule resid {h2_molecule} not found in {replica_name}. Skipping this H2 molecule.")
            continue

        # Create a selection that includes both the protein and the specified H2 molecule
        selection = protein + h2_molecule_atoms

        # Set up the reference Universe (first iteration)
        if ref_universe is None:
            # Load the .gro file of the first replica as the reference
            ref_universe = mda.Universe(top_path)  # This only loads the topology (gro file)
            ref_atoms = ref_universe.select_atoms("backbone")  # Select the backbone atoms for alignment
            ref_universe.atoms.positions  # Ensure the atoms positions are loaded

            # Save the reference PDB (optional)
            ref_universe.atoms.write(os.path.join(output_dir, "reference.pdb"))

        # Apply alignment using the reference's backbone
        alignment = align.AlignTraj(u, ref_universe, select="backbone", in_memory=True)
        alignment.run()  # Perform the alignment

        # Write the aligned trajectory
        aligned_traj_path = os.path.join(output_dir, f'{replica_name}_H2_{h2_molecule}.xtc')
        with mda.Writer(aligned_traj_path, selection.n_atoms) as writer:
            for ts in u.trajectory:
                writer.write(selection)  # Write the aligned frame

print("Processing and alignment complete for all replicas.")

