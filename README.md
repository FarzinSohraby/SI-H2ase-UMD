# SI-H2ase-UMD
This is the Supporting Information for the paper of UMD simulations of H2ase

1-Distance-H2-MNI-allreps.py
	This code has been used for calculating the distance values of each of the 100 H2 molecules in all replicas. A threshold of 5 angstrom was then used to find the binding events.
	
2- pbc-align-xtc.py
	This code has been used for aligning the 75 replicas to one reference (the gro file).
	
3- Df-H2-10ps-FEL.ipynb
	This notebook was used for the dimentionality reduction with TICA and then to plot the Free Energy Landscape figures for Df H2ase.
	
4- Mdg-H2-10ps-FEL.ipynb
	This notebook was used for the dimentionality reduction with TICA and then to plot the Free Energy Landscape figures for Mdg H2ase.
	
5-FPT.py
	This code has been used for calculating the First passsage times based on the 5 angstrom threshold.
	
6-Chi2-test-Df.ipynb
	This notebook was used for the chi-square statistical test for the binding and unbinding pathways of Df H2ase.

7-Chi2-test-Mdg.ipynb
	This notebook was used for the chi-square statistical test for the binding and unbinding pathways of Mdg H2ase.

8 and 9- Pymol session files containing the Df and Mdg H2ase tunnels indentified by CAVER plug-in tool in pymol.
