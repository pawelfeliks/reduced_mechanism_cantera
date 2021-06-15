import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.yaml')
initial_state = 1200, 5 * ct.one_atm, 'CH4:0.35, O2:1.0, N2:3.76'

# Run a simulation with the full mechanism
gas.TPX = initial_state
r = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([r])

tt = []
TT = []
t = 0.0
# Rmax is the maximum relative reaction rate at any timestep
Rmax = np.zeros(gas.n_reactions)

print(gas.n_reactions)

while t < 0.02:
    t = sim.step()
    tt.append(1000 * t)
    TT.append(r.T)
    rnet = abs(gas.net_rates_of_progress)
    rnet /= max(rnet)
    Rmax = np.maximum(Rmax, rnet)

plt.plot(tt, TT, label='K=53, R=325', color='k', lw=3, zorder=100)

# Get the reaction objects, and sort them so the most active reactions are first
R = sorted(zip(Rmax, gas.reactions()), key=lambda x: -x[0])

#print('\n NEW ORDER:')

#print(gas.species())
#print(gas.reactions())
#C = plt.cm.winter(np.linspace(0, 1, 1))
N = 60  #number of avtive reactions
for i in range(N):
    # Get the 40 most active reactions
    reactions = [r[1] for r in R[:N]]
    
    # find the species involved in these reactions. At a minimum, include all
    # species in the reactant mixture
    species_names = {'N2', 'CH4', 'O2'}
    for reaction in reactions:
        species_names.update(reaction.reactants)
        species_names.update(reaction.products)
    
    # Get the species objects
    species = [gas.species(name) for name in species_names]
    # create the new reduced mechanism
    gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                       species=species, reactions=reactions)
    
    gas2.TPX = initial_state
    r = ct.IdealGasConstPressureReactor(gas2)
    sim = ct.ReactorNet([r])
    gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                       species=species, reactions=reactions)
     
    # Re-run the ignition problem with the reduced mechanism
    gas2.TPX = initial_state
    r = ct.IdealGasConstPressureReactor(gas2)
    sim = ct.ReactorNet([r])

    t = 0.0

    tt = []
    TT = []
    while t < 0.02:
        t = sim.step()
        tt.append(1000 * t)
        TT.append(r.T)

  
plt.plot(tt, TT, label='K={0}, R={1}'.format(gas2.n_species, N))
plt.xlabel('Time (ms)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper left')
plt.title('Reduced mechanism ignition delay times\n'
              'K: number of species; R: number of reactions')
plt.xlim(0, 20)
plt.tight_layout()



plt.show()

'''
print(gas.species())
print('\n REDUCED SPECIES:')
print(gas2.species())
print('\n REDUCED Reactions:')
print(gas2.reactions())
'''
