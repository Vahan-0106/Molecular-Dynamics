from openmm import unit
from openmm.app import *
from openmm import *
from sys import stdout
import mdtraj as md
import os, sys

pdb = PDBFile("aki_centered.pdb")

forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

modeller = Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
modeller.addHydrogens(forcefield)
modeller.addSolvent(
    forcefield=forcefield,
    padding=1.0 * unit.nanometer,
    ionicStrength=0.15 * unit.molar,
    neutralize=True,
)

with open("solvated_system.pdb", "w") as fh:
    PDBFile.writeFile(modeller.topology, modeller.positions, fh)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=HBonds,
)

integrator = LangevinMiddleIntegrator(
    10 * unit.kelvin,           
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds,
)

system.addForce(MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin, 25))

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy()

dcd = DCDReporter("trajectory.dcd", 200)        
simulation.reporters.append(dcd)
simulation.reporters.append(
    StateDataReporter(
        stdout, 200, step=True, temperature=True, totalEnergy=True, density=True
    )
)

for T in range(15, 305, 5):
    print(f"Heating to {T} K …")
    integrator.setTemperature(T * unit.kelvin)
    simulation.step(100)

simulation.step(30000)
dcd._out.close()                                     

state = simulation.context.getState(getPositions=True)
with open("final_boxed_structure.pdb", "w") as fh:
    PDBFile.writeFile(simulation.topology, state.getPositions(), fh)

traj = md.load("trajectory.dcd", top="solvated_system.pdb")
traj.save_xtc("trajectory.xtc")

print("Done.   DCD size:", os.path.getsize("trajectory.dcd"), "bytes")
