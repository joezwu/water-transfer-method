import MDAnalysis
from MDAnalysis.analysis import distances
import sys

basename = sys.argv[1]
cutoff = float(sys.argv[2])
ligand_hselection = sys.argv[3]
ligand_hselection2 = sys.argv[4]

traj = basename + "_md.xtc"

u = MDAnalysis.Universe(basename + "_md.pdb", traj)

ow = u.select_atoms('( resname HOH and name O) and around %f ( %s )' % (cutoff, ligand_hselection), periodic = True, updating=True)
#ow = u.select_atoms('( resname HOH and name O) and (around %f ( %s ) or around %f ( %s ))' % (cutoff, ligand_hselection, cutoff, ligand_hselection2), periodic = True, updating=True)
#ow = u.select_atoms('(( resname HOH and name O) and (around %f ( %s ))) or (( resname HOH and name O) and (around %f ( %s )))' % (cutoff, ligand_hselection, cutoff, ligand_hselection2), periodic = True, updating=True)

#ow = u.select_atoms('( resname HOH and name O) and around 4.0 (resid 265 and resname LIG and name N2)', periodic = True, updating=True)

for ts in u.trajectory:     # iterate through all frames
    #print("Num Waters:", len(ow))
    print(" %d " % len(ow))
