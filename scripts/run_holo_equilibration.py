import sys
import openmm as mm
import openmm.app as app
import openmm.unit as unit

print("[*] Cargando topología ternaria SIRT6_c1 + Zn + MDL800...")
prmtop = app.AmberPrmtopFile('results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop')
inpcrd = app.AmberInpcrdFile('results/holo_mdl800/SIRT6_Holo_MDL800_solv.inpcrd')

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=app.HBonds,
    rigidWater=True
)

# Escudo armónico sobre Esqueleto + Clúster Cuántico Zn-S4 + Átomos pesados de MDL800
restraint = mm.CustomExternalForce("0.5*k_rest*periodicdistance(x, y, z, x0, y0, z0)^2")
restraint.addGlobalParameter("k_rest", 10.0*unit.kilocalories_per_mole/unit.angstroms**2)
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

restrained_count = 0
for atom in prmtop.topology.atoms():
    resname = atom.residue.name
    is_backbone = atom.name in ('CA', 'C', 'N')
    is_zn_cluster = resname in ('CM1', 'CM2', 'CM3', 'CM4', 'ZN1') and (atom.element is None or atom.element.symbol != 'H')
    is_ligand_heavy = resname == 'MDL' and (atom.element is None or atom.element.symbol != 'H')
    if is_backbone or is_zn_cluster or is_ligand_heavy:
        pos = inpcrd.positions[atom.index]
        restraint.addParticle(atom.index, [pos[0], pos[1], pos[2]])
        restrained_count += 1

system.addForce(restraint)
print(f"[*] Escudo armónico activado sobre {restrained_count} átomos clave (Backbone + Zn-S4 + MDL800).")

integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}

simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print("[*] Minimizando choques estéricos de solvatación (Max 2000 iteraciones)...")
simulation.minimizeEnergy(maxIterations=2000)

simulation.context.setVelocitiesToTemperature(300*unit.kelvin, 42)
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, 25000, step=True, potentialEnergy=True, temperature=True, speed=True
))
simulation.reporters.append(app.DCDReporter('results/holo_mdl800/SIRT6_Holo_equil.dcd', 25000))

print("[*] Fase 1: Termalización NVT a 300 K (100 ps, k = 10.0 kcal/mol/Å²)...")
simulation.step(50000)

print("[*] Fase 2: Activando Barostato de Monte Carlo (1 bar) y relajando escudo gradualmente (1.0 ns)...")
system.addForce(mm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
simulation.context.reinitialize(preserveState=True)

for k_val in [5.0, 2.5, 1.0, 0.5, 0.0]:
    print(f"    -> Relajación NPT con k_rest = {k_val} kcal/mol/Å² (200 ps)...")
    simulation.context.setParameter(
        "k_rest", k_val*unit.kilocalories_per_mole/unit.angstroms**2
    )
    simulation.step(100000)

simulation.saveState('results/holo_mdl800/SIRT6_Holo_equil.xml')
print("[✓] Equilibración Holo-MD completada exitosamente.")
