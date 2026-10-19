import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys
import os

print("[*] Cargando topología AMBER...")
prmtop = app.AmberPrmtopFile('results/mcpb/SIRT6_Zn_solv.prmtop')

system = prmtop.createSystem(nonbondedMethod=app.PME, 
                             nonbondedCutoff=1.0*unit.nanometers,
                             constraints=app.HBonds, 
                             rigidWater=True)
system.addForce(mm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))

integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}

simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)

chk_file = 'results/mcpb/SIRT6_Zn_prod.chk'
total_steps = 25000000

if os.path.exists(chk_file):
    print("[*] Resurrección térmica: Reanudando desde el último Checkpoint...")
    simulation.loadCheckpoint(chk_file)
    current_step = simulation.currentStep
else:
    print("[*] Extrayendo física pura del estado equilibrado...")
    with open('results/mcpb/SIRT6_Zn_equil.xml', 'r') as f:
        state = mm.XmlSerializer.deserialize(f.read())
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())
    simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    current_step = 0
    # Eliminamos restos de DCD corrupto si estamos empezando de cero
    if os.path.exists('results/mcpb/SIRT6_Zn_prod_50ns.dcd'):
        os.remove('results/mcpb/SIRT6_Zn_prod_50ns.dcd')

steps_left = total_steps - current_step
append_mode = current_step > 0

simulation.reporters.append(app.StateDataReporter(sys.stdout, 25000, step=True, 
    potentialEnergy=True, temperature=True, speed=True))
simulation.reporters.append(app.DCDReporter('results/mcpb/SIRT6_Zn_prod_50ns.dcd', 25000, append=append_mode))
# Checkpoint cada 250,000 pasos (0.5 ns). Si falla, perdemos máximo 10 minutos.
simulation.reporters.append(app.CheckpointReporter(chk_file, 250000))

print(f"[*] Detonando Dinámica Molecular. Pasos restantes: {steps_left}...")
if steps_left > 0:
    simulation.step(steps_left)
    simulation.saveState('results/mcpb/SIRT6_Zn_prod_final.xml')
    print("[✓] Muestreo conformacional completado.")
else:
    print("[✓] La trayectoria ya está completa.")
