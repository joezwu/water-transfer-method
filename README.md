# ATM Water-Transfer-Method
Simulation input files for the paper "Ligand-Water Alchemical Transfer to Model Bridging Solvent in Relative Binding Free Energy Calculations" by Joe Z Wu and Emilio Gallicchio. Pre-print can be accessed at https://arxiv.org/abs/.

The water-transfer simulations require AToM-OpenMM v3.5.0.

Water-Transfer RBFE allows for the relative binding estimation between two ligands that bind with a variable number of bridging waters.

Included in this repository are simulation input files that should allow reproduction of the RBFE estimates obtained in the above publication. Three simulations were conducted: 
(1) Water-Transfer,
(2) regular/conventional RBFE starting with the bridging water, and 
(3) regular/conventional RBFE starting without the bridging water. 

For the Water-Transfer simulation, an additional water-counting simulation was performed on the ligand in solution to correct for adding the Water-Transfer tether.
