import os
import sys
import struct
import openmm as mm
import openmm.app as app
import openmm.unit as unit

print("[*] Cargando topología AMBER Holo (SIRT6_c1 + Zn + MDL800)...")
prmtop = app.AmberPrmtopFile('results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop')

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=app.HBonds,
    rigidWater=True
)
system.addForce(mm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))

integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}

simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)

chk_file = 'results/holo_mdl800/SIRT6_Holo_prod.chk'
dcd_file = 'results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd'
dcd_bak = 'results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd.bak'
total_steps = 25000000
report_interval = 25000

if os.path.exists(chk_file) and (os.path.exists(dcd_file) or os.path.exists(dcd_bak)):
    if not os.path.exists(dcd_file):
        os.link(dcd_bak, dcd_file)

    print("[*] Resurrección térmica: Reanudando Holo-MD desde el último Checkpoint...")
    simulation.loadCheckpoint(chk_file)
    current_step = simulation.context.getState().getStepCount()
    simulation.currentStep = current_step

    # En OpenMM DCDFile, la cabecera mide exactamente 276 bytes y NSET reside en el byte 8
    expected_frames = current_step // report_interval
    header_total = 276
    frame_bytes = 80 + 12 * prmtop.topology.getNumAtoms()
    with open(dcd_file, 'r+b') as f:
        f.seek(8)
        f.write(struct.pack('<i', expected_frames))
        f.truncate(header_total + expected_frames * frame_bytes)
    print(f"[✓] Estado sincronizado en el paso {current_step} ({expected_frames} fotogramas DCD verificados).")
else:
    print("[*] Extrayendo física pura del estado Holo equilibrado (desde paso 0)...")
    with open('results/holo_mdl800/SIRT6_Holo_equil.xml', 'r') as f:
        state = mm.XmlSerializer.deserialize(f.read())
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())
    simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    current_step = 0
    if os.path.exists(dcd_file):
        os.remove(dcd_file)

steps_left = total_steps - current_step
append_mode = current_step > 0

simulation.reporters.append(app.StateDataReporter(
    sys.stdout, report_interval, step=True, potentialEnergy=True, temperature=True, speed=True
))
simulation.reporters.append(app.DCDReporter(dcd_file, report_interval, append=append_mode))
simulation.reporters.append(app.CheckpointReporter(chk_file, 250000))

print(f"[*] Ejecutando Dinámica Molecular Holo. Pasos restantes: {steps_left} ({steps_left*0.002/1000:.1f} ns)...")
if steps_left > 0:
    simulation.step(steps_left)
    simulation.saveState('results/holo_mdl800/SIRT6_Holo_prod_final.xml')
    if os.path.exists(dcd_bak):
        os.remove(dcd_bak)
    print("[✓] Muestreo conformacional Holo-MD (50 ns) completado.")
