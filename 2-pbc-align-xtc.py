import os
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align

# Paths
csv_path = '/media/obiwan/Elements1/NewHyd-WT-100H2-R2/summary_distances.csv'
replica_dir = '/media/obiwan/Elements1/NewHyd-WT-100H2-R2'
output_dir = '/media/obiwan/Elements1/NewHyd-WT-100H2-R2/Aligned-xtc-all-H2'

# Load CSV data (used for matching replica names, even though ligands are not selected explicitly)
csv_data = pd.read_csv(csv_path)

# Set up directories for saving processed trajectories
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Initialize reference Universe (first replica)
ref_universe = None
ref_atoms = None

# Iterate over each replica and process trajectories
for i in range(24, 30):
    replica_name = f'Replica{i}'  # Match exact format in the CSV
    traj_path = os.path.join(replica_dir, replica_name, 'NewHyd-WT-H2-Bind-250ns-pbc-ligs.xtc')  # Trajectory path
    top_path = '/media/obiwan/Elements1/NewHyd-WT-100H2-R2/NewHyd-WT-H2-Bind-250ns-ligs.gro'  # Topology file

    # Check if the replica exists in the CSV
    if not csv_data[csv_data['Replica'] == replica_name].empty:
        print(f"Processing {replica_name}...")
    else:
        print(f"{replica_name} not found in CSV data. Skipping...")
        continue

    # Load the trajectory
    try:
        u = mda.Universe(top_path, traj_path)
    except Exception as e:
        print(f"Error loading trajectory or topology for {replica_name}: {e}")
        continue

    # Select the protein and all ligands
    protein = u.select_atoms('protein')
    ligands = u.select_atoms('not protein')  # Select everything that is not protein
    selection = protein + ligands  # Combine the protein and ligand selections

    # Set up the reference Universe (first iteration)
    if i == 24:
        # Load the .gro file of the first replica as the reference
        ref_universe = mda.Universe(top_path)  # Load topology
        ref_atoms = ref_universe.select_atoms("backbone")  # Select the backbone atoms for alignment
        ref_universe.atoms.positions  # Ensure the positions are loaded

        # Save the reference PDB (optional)
        ref_universe.atoms.write(os.path.join(output_dir, "reference.pdb"))

    # Apply alignment using the reference's backbone
    alignment = align.AlignTraj(u, ref_universe, select="backbone", in_memory=True)
    alignment.run()  # Perform the alignment

    # Write the aligned trajectory
    aligned_traj_path = os.path.join(output_dir, f'{replica_name}.xtc')
    with mda.Writer(aligned_traj_path, selection.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(selection)  # Write the aligned frame

print("Processing and alignment complete for all replicas.")

