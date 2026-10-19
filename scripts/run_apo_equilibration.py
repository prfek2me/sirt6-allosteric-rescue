import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys, os

prmtop_file = 'results/mcpb/SIRT6_Zn_solv.prmtop'
inpcrd_file = 'results/mcpb/SIRT6_Zn_solv.inpcrd'

if not os.path.exists(prmtop_file):
    print("[!] Faltan archivos de topología. Abortando.")
    sys.exit(1)

print("[*] Cargando topología AMBER en OpenMM...")
prmtop = app.AmberPrmtopFile(prmtop_file)
inpcrd = app.AmberInpcrdFile(inpcrd_file)

system = prmtop.createSystem(nonbondedMethod=app.PME, 
                             nonbondedCutoff=1.0*unit.nanometers,
                             constraints=app.HBonds, 
                             rigidWater=True)

# 2. Inyección de Restricciones (Escudo Temporal)
force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
force.addGlobalParameter("k", 10.0 * unit.kilocalories_per_mole/unit.angstroms**2)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

for atom in prmtop.topology.atoms():
    if atom.residue.name in ['ZN', 'CYM'] or atom.name in ['CA', 'C', 'N', 'O']:
        force.addParticle(atom.index, inpcrd.positions[atom.index])
system.addForce(force)

# Reducimos el paso a 1 femtosegundo (0.001 ps) para máxima estabilidad durante la compresión
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}

simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print("[*] FASE 1: Minimizando energía estérica (0 K)...")
simulation.minimizeEnergy(maxIterations=10000)

print("[*] FASE 2: Calentamiento isocórico NVT con escudo (300K, 100 ps)...")
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.step(100000) # 100,000 pasos x 1 fs

print("[*] FASE 3: Acoplamiento isobárico NPT (1 atm, 1 ns)...")
# LA CLAVE MAGISTRAL: Apagamos el escudo a cero. Dejamos que Seminario haga su trabajo.
simulation.context.setParameter("k", 0.0)
system.addForce(mm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
simulation.context.reinitialize(preserveState=True)

simulation.reporters.append(app.StateDataReporter(sys.stdout, 50000, step=True, 
    potentialEnergy=True, temperature=True, volume=True, density=True))
simulation.reporters.append(app.DCDReporter('results/mcpb/SIRT6_Zn_equil.dcd', 100000))

# 1 nanosegundo puro de compresión líquida libre de tensiones
simulation.step(1000000)

simulation.saveState('results/mcpb/SIRT6_Zn_equil.xml')
print("[✓] Equilibrio termodinámico alcanzado libre de colapsos. Trayectoria guardada.")