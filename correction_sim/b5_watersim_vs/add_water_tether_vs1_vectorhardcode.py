# Usage: python add_water_tether_vs1_vectorhardcode.py --systemXMLinFile b5_sys_noVS.xml --systemXMLoutFile b5_sys.xml --receptorPDBinFile b5_noVS.pdb --systemPDBoutFile b5.pdb --virtualsite "21 22 23"
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
parser.add_argument('--virtualsite',  required=True,  type=str, default=None,
        help='Adds Three-Particle-Average Site: "index1 index2 index3"')

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
vs1 = args['virtualsite']

#adds topology information for virtual site
chaincomplex1 = pdbrcpt.topology.addChain("G")
residuecomplex1 = pdbrcpt.topology.addResidue("VS1", chaincomplex1)
element=Element.getByAtomicNumber(1)


#calculates virtual site vector from ThreeParticleAverageSite
#indexes are from left to right, with the middle index being the center of the three points
def addVS1(virtualsite):
    vs_indexes = [int(virtualsite.split()[0]), int(virtualsite.split()[1]), int(virtualsite.split()[2])]

    r0_v = np.asarray(rcpt_positions[vs_indexes[0]]/unit.nanometer) 
    r1_v = np.asarray(rcpt_positions[vs_indexes[1]]/unit.nanometer) 
    r2_v = np.asarray(rcpt_positions[vs_indexes[2]]/unit.nanometer)
    #d = np.linalg.norm(np.asarray(rcpt_positions[int(watom1)]/unit.nanometer)-r1_v)

    #h = d / np.sqrt(4*np.dot(r1_v, r1_v) - 4*np.dot(r1_v, r0_v + r2_v) + np.dot(r0_v + r2_v, r0_v + r2_v))

    #w0 = -1*h
    #w1 = 1+2*h
    #w2 = -1*h
    w0 = -2.3145269711491294
    w1 = 5.629053942298259
    w2 = -2.3145269711491294
    rv = (w0*r0_v + w1*r1_v + w2*r2_v) * unit.nanometer / unit.angstrom
    rv_vec3 = Vec3(rv[0], rv[1], rv[2]) * angstrom

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

with open(args['systemXMLoutFile'], 'w') as output:
    output.write(XmlSerializer.serialize(system))


if pdboutfile is not None:
    PDBFile.writeFile(pdbrcpt.topology, pdbrcpt.positions,
                      open(pdboutfile,'w'))
