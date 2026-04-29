from __future__ import print_function

import openmm as mm
from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

#the multiple-time step integrator does not have a setTemperature() method
def setTemperature(self, temperature):
    self.setGlobalVariableByName('kT', MOLAR_GAS_CONSTANT_R*temperature)
MTSLangevinIntegrator.setTemperature = setTemperature

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "c3d"

#load system
pdb = PDBFile(jobname + ".pdb")
topology = pdb.topology
positions = pdb.positions
boxVectors = topology.getPeriodicBoxVectors()
systemfile = jobname + "_sys.xml"
#HMASS is set in the system file
with open(systemfile) as input:
    system = XmlSerializer.deserialize(input.read())

#Set up Langevin integrator
initial_temperature = 50 * kelvin
final_temperature = 300 * kelvin
temperature = initial_temperature
frictionCoeff = 0.05 / picosecond
MDstepsize = 0.001 * picosecond
barostat = MonteCarloBarostat(1*bar, final_temperature)
saved_barostat_frequency = barostat.getFrequency()
barostat.setFrequency(0)#disabled
system.addForce(barostat)
integrator = LangevinIntegrator(temperature, frictionCoeff, MDstepsize)
integrator.setConstraintTolerance(0.00001)

#platform_name = 'OpenCL'
platform_name = 'CUDA'
#platform_name = 'HIP'
platform = Platform.getPlatformByName(platform_name)
properties = {}
properties["Precision"] = "mixed"

simulation = Simulation(topology, system, integrator, platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
print(boxVectors)
simulation.context.setPositions(positions)
if boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*boxVectors)

print("Potential energy before minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Energy minimizing the system ...")
simulation.minimizeEnergy()

print("Potential energy after minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

#saves checkpoint
simulation.saveState(jobname + '_min.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_min.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

print("Thermalization ...")

totalSteps = 150000
steps_per_cycle = 5000
number_of_cycles = int(totalSteps/steps_per_cycle)
delta_temperature = (final_temperature - initial_temperature)/number_of_cycles
simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True))

#MD with temperature ramp
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)
    #prepare system for new temperature
    temperature = temperature + delta_temperature
    integrator.setTemperature(temperature)

#saves checkpoint
simulation.saveState(jobname + '_therm.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_therm.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
    
print("NPT equilibration ...")
barostat.setFrequency(saved_barostat_frequency)#enabled

#MD at constant pressure
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)

#saves checkpoint
simulation.saveState(jobname + '_npt.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_npt.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

print("NVT production ...")
barostat.setFrequency(0)#disabled

#30,000,000 fs is equal to 30000000 steps * (1fs/step)
prod_totalSteps = 30000000
prod_steps_per_cycle = 5000

simulation.reporters.append(XTCReporter(jobname + "_md.xtc", prod_steps_per_cycle ))

#MD at constant volume
simulation.step(prod_totalSteps)

#saves checkpoint
simulation.saveState(jobname + '_md.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_md.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
    
end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
