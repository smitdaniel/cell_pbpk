from typing import Callable

from matplotlib import pyplot as plt

from src.cell import Cell
from src.compartment import Compartment, SpecieType
from src.pathway import SimplePathway

if __name__ == "__main__":

    """
    Parameters
    """

    # how many molecules of each specie are injected on cmp.inject() call
    SPECIE_INJECTION = {
        SpecieType.A: 1000,
        SpecieType.cA: 1000,
        SpecieType.B: 500,
        SpecieType.cB: 75,
        SpecieType.C: 90
    }

    # receptor occupancy ratio to switch to ON/True state
    # the pathway is a logical model:
    # (A > 0.25)    --| A and not B  --|
    # (B > 0.4)     --|                |-- AB or not BC
    # (C > 0.1)     --| B or not C   --|
    PATHWAY = SimplePathway({
        SpecieType.A: 0.25,
        SpecieType.B: 0.4,
        SpecieType.C: 0.1})
    # number of receptors on the cell
    RECEPTOR_COUNT = {
        SpecieType.A: 1000,
        SpecieType.B: 500,
        SpecieType.C: 5000
    }

    # time clicks
    TICKTOCKS: int = 150000
    do_inject: Callable[[int], bool] = lambda t: t % 10000 == 0

    """
    Objects
    """

    # kinetic rates and decay of species are defined in compartment.py
    cmp = Compartment(SPECIE_INJECTION)
    cmp.inject()

    # cell with a given PATHWAY and an output rate per time step
    cell = Cell(PATHWAY, 0.5, RECEPTOR_COUNT)

    # timesteps containers
    time = []
    occupancy_ratios = {st: [] for st in SpecieType}
    compartment_concentrations = {st: [] for st in SpecieType}

    for t in range(TICKTOCKS):
        time.append(t)
        if do_inject(t):
            cmp.inject()
        # compute species decay in the compartment
        cmp.update_species_decay()
        # compute on/off receptor binding changes on the cell
        # release/absorb species from the compartment
        cell.update_occupancy(comp=cmp)
        # depending on the pathway, produce SpecieType.B, release into compartment
        cmp.species[SpecieType.B] += cell.production
        # record occupancy
        for st in occupancy_ratios.keys():
            occupancy_ratios[st].append(cell.get_occupancy_ratio(st))
        # record concentration of SpecieType.B in the compartment
        compartment_concentrations[SpecieType.B].append(cmp.get_concentration(SpecieType.B))

    # plot
    fig, ax = plt.subplots(2)
    clr = {SpecieType.A: 'r', SpecieType.B: 'b', SpecieType.C: 'g',
           SpecieType.cB: 'y', SpecieType.cA: 'k'}
    # plot thresh lines
    ax[0].axhline(0.25, label='thresh A', c='r', linestyle='--')
    ax[0].axhline(0.4, label='thresh B', c='b', linestyle='--')
    ax[0].axhline(0.1, label='thresh C', c='g', linestyle='--')
    ax[0].set_title('Occupancy fraction')
    ax[1].set_title('SpeciesB concentration in the Compartment')
    # plot occupancy ratios
    for st, aor in occupancy_ratios.items():
        if st == SpecieType.cA:
            continue
        ax[0].plot(time, aor, label=str(st), c=clr[st])
    # plot compartment concentration
    ax[1].plot(time, compartment_concentrations[SpecieType.B], label=f"conc-{SpecieType.B}")
    ax[0].legend()
    plt.show()
