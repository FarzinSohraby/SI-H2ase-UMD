Supporting Information for the Paper on UMD Simulations of H₂ase

Distance-H2-MNI-allreps.py:
This code was used to calculate the distance values of each of the 100 H₂ molecules in all replicas. A threshold of 5 Ångström was then used to identify binding events.

pbc-align-xtc.py:
This code was used to align the 75 replicas to a single reference (the .gro file).

FPT.py:
This code was used to calculate the First Passage Times based on the 5 Å threshold.

Chi2-test-Df.ipynb:
This notebook was used for the chi-square statistical test for the binding and unbinding pathways of Df H₂ase.

Chi2-test-Mdg.ipynb:
This notebook was used for the chi-square statistical test for the binding and unbinding pathways of Mdg H₂ase.

Pymol session files:
These contain the Df and Mdg H₂ase tunnels identified using the CAVER plug-in tool in PyMOL.

MSM-model-Building-Df.ipynb:
This notebook contains the codes used to build and validate the MSM and also compute the rates for Df hydrogenase.

MSM-model-Building-Mdg.ipynb:
This notebook contains the codes used to build and validate the MSM and also compute the rates for Mdg hydrogenase.
