# Usage: python add_water_tether_vs2.py --systemXMLinFile btk-s8-s11_sys_noVS.xml --systemXMLoutFile btk-s8-s11_sys.xml --receptorPDBinFile btk-s8-s11_noVS.pdb --systemPDBoutFile btk-s8-s11.pdb --virtualsite "5404 5414 5405" --waterOindex1 "4849" --waterOindex2 ""
# Emilio Gallicchio, 5/2023 adapted from code by Bill Swope, 11/2021

############################################
#                                          #
#   IMPORTS                                #
#                                          #
############################################

import os, sys
import numpy as np
from datetime import datetime
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm import XmlSerializer 
from openmm import Vec3
from openmm.app import PDBReporter, StateDataReporter, PDBFile
from openmm.app import ForceField, Modeller
from openmm.app import PME, HBonds

from openmm.unit import Quantity
from openmm.unit import angstrom, nanometer, nanometers, picoseconds, amu, kilocalorie_per_mole, kilojoule_per_mole
from openmm import unit
from openmm import CustomBondForce
from sys import stdout

#from pdbfixer import PDBFixer

# OpenFF components for SystemGenerator (for input entire ssytem in pdb file)
from openmm import app
from openmmforcefields.generators import SystemGenerator
from openmm import ThreeParticleAverageSite

import re
from openmm.app import Element

#####################################################
#                                                   #
#   PASS COMMAND LINE ARGUMENTS TO LOCAL VARIABLES  #
#                                                   #
#####################################################

whatItDoes = ""

program_start_timer = time()
parser = argparse.ArgumentParser(description=whatItDoes)

# Required input
parser.add_argument('--receptorPDBinFile', required=True,  type=str, default=None,
                    help='Receptor PDB format (.pdb)')
parser.add_argument('--systemPDBoutFile', required=True,  type=str, default=None,
                    help='Complex PDB format (.pdb)')
parser.add_argument('--systemXMLinFile',  required=True,  type=str, default=None,
                    help='Input XML file where to save the System')
parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None,
                    help='Output XML file where to save the System')
parser.add_argument('--virtualsite1',  required=False,  type=str, default=None,
        help='Adds Three-Particle-Average Site: "index1, index2, index3"')
parser.add_argument('--virtualsite2',  required=False,  type=str, default=None,
        help='Adds second Three-Particle-Average Site: "index1, index2, index3"')
parser.add_argument('--waterOindex1',  required=True,  type=str, default=None,
        help='Adds index of water oxygen: "indexO_1"')
parser.add_argument('--waterOindex2',  required=True,  type=str, default=None,
        help='Adds index of water oxygen: "indexO_2"')


# Arguments that are flags
parser.add_argument('--verbose', required=False, action='store_true',
                    help='Get more output with this flag')

args = vars(parser.parse_args())

with open(args['systemXMLinFile']) as input:
    system = XmlSerializer.deserialize(input.read())

pdboutfile = args['systemPDBoutFile']

receptorfile = args['receptorPDBinFile']
pdbrcpt = PDBFile(receptorfile)
rcpt_positions = pdbrcpt.positions

#complexfile = args['systemPDBinFile']
#pdbcomplex = PDBFile(complexfile)

#three point VMD indexes of ligand atoms to form virtual site
vs1 = args['virtualsite1']
vs2 = args['virtualsite2']
#index of oxygen of first tethered water (VMD index)
watom1 = args['waterOindex1']
#index of oxygen of second tethered water (VMD index)
watom2 = args['waterOindex2']

#adds topology information for virtual site
chaincomplex1 = pdbrcpt.topology.addChain("H")
residuecomplex1 = pdbrcpt.topology.addResidue("VS1", chaincomplex1)
chaincomplex2 = pdbrcpt.topology.addChain("I")
residuecomplex2 = pdbrcpt.topology.addResidue("VS2", chaincomplex2)
element=Element.getByAtomicNumber(1)


#calculates virtual site vector from ThreeParticleAverageSite
def addVS1(virtualsite):
    vs_indexes = [int(virtualsite.split()[0]), int(virtualsite.split()[1]), int(virtualsite.split()[2])]
    print(vs_indexes)
    r0_v = np.asarray(rcpt_positions[vs_indexes[0]]/unit.nanometer) 
    r1_v = np.asarray(rcpt_positions[vs_indexes[1]]/unit.nanometer) 
    r2_v = np.asarray(rcpt_positions[vs_indexes[2]]/unit.nanometer)
    d = np.linalg.norm(np.asarray(rcpt_positions[int(watom1)]/unit.nanometer)-r1_v)
    print(int(watom1)) 
    h = d / np.sqrt(4*np.dot(r1_v, r1_v) - 4*np.dot(r1_v, r0_v + r2_v) + np.dot(r0_v + r2_v, r0_v + r2_v))
    w0 = -1*h
    w1 = 1+2*h
    w2 = -1*h
    rv = (w0*r0_v + w1*r1_v + w2*r2_v) * unit.nanometer / unit.angstrom
    rv_vec3 = Vec3(rv[0], rv[1], rv[2]) * angstrom
    print(rv_vec3)
    p = system.addParticle(0.0)
    system.setVirtualSite(p, ThreeParticleAverageSite(int(virtualsite.split()[0]), int(virtualsite.split()[1]), int(virtualsite.split()[2]), w0, w1, w2))
    pdbrcpt.positions.append(rv_vec3)

nbpattern = re.compile(".*Nonbonded.*")
for i in range(system.getNumForces()):
    if nbpattern.match(str(type(system.getForce(i)))):
        nbforce = system.getForce(i)
        break
#add virtual particle to nbforce
nbforce.addParticle(0.0, 0.1, 0.0)

addVS1(vs1)
vs1_atom = pdbrcpt.topology.addAtom("VS1", element, residuecomplex1)

#calculates virtual site vector from ThreeParticleAverageSite
def addVS2(virtualsite):
    vs_indexes = [int(virtualsite.split()[0]), int(virtualsite.split()[1]), int(virtualsite.split()[2])]
    print(vs_indexes)
    r0_v = np.asarray(rcpt_positions[vs_indexes[0]]/unit.nanometer) 
    r1_v = np.asarray(rcpt_positions[vs_indexes[1]]/unit.nanometer)
    r2_v = np.asarray(rcpt_positions[vs_indexes[2]]/unit.nanometer)
    d = np.linalg.norm(np.asarray(rcpt_positions[int(watom2)]/unit.nanometer)-r1_v)
    print(int(watom2)) 
    h = d / np.sqrt(4*np.dot(r1_v, r1_v) - 4*np.dot(r1_v, r0_v + r2_v) + np.dot(r0_v + r2_v, r0_v + r2_v))
    w0 = -1*h
    w1 = 1+2*h
    w2 = -1*h
    rv = (w0*r0_v + w1*r1_v + w2*r2_v) * unit.nanometer / unit.angstrom
    rv_vec3 = Vec3(rv[0], rv[1], rv[2]) * angstrom
    print(rv_vec3)
    r = system.addParticle(0.0)
    system.setVirtualSite(r, ThreeParticleAverageSite(int(virtualsite.split()[0]), int(virtualsite.split()[1]), int(virtualsite.split()[2]), w0, w1, w2))
    pdbrcpt.positions.append(rv_vec3)
    
nbpattern = re.compile(".*Nonbonded.*")
for i in range(system.getNumForces()):
    if nbpattern.match(str(type(system.getForce(i)))):
        nbforce = system.getForce(i)
        break
#add virtual particle to nbforce
nbforce.addParticle(0.0, 0.1, 0.0)

addVS2(vs2)
vs2_atom = pdbrcpt.topology.addAtom("VS2", element, residuecomplex2)

#water tethering potential
wbforce = CustomBondForce("0.5*kfw*max(0,r-tolw)^2")
wbforce.addPerBondParameter("kfw")
wbforce.addPerBondParameter("tolw")
#index of oxygen of tethered water (VMD index)
watom1 = int(args['waterOindex1'])
watom2 = int(args['waterOindex2'])
#index of ligand atom tethered to water (VMD index)
lig_watom1 = vs1_atom.index
lig_watom2 = vs2_atom.index
#sets the force constant and the tolerance parameters of the tethering potential
kfw = ( 25. * kilocalorie_per_mole/angstrom**2 )/(kilojoule_per_mole/nanometer**2)
tolw = ( 3.0 * angstrom )/nanometer
#adds the tether to the system
wbforce.addBond(watom1, lig_watom1, [kfw, tolw])
wbforce.addBond(watom2, lig_watom2, [kfw, tolw])
system.addForce(wbforce)

with open(args['systemXMLoutFile'], 'w') as output:
    output.write(XmlSerializer.serialize(system))


if pdboutfile is not None:
    PDBFile.writeFile(pdbrcpt.topology, pdbrcpt.positions,
                      open(pdboutfile,'w'))
