from dataclasses import dataclass
from functools import cached_property
from typing import Dict

from src.compartment import SpecieType, Compartment
from src.pathway import SimplePathway


@dataclass
class Occupancy:
    free: float  # free receptors
    ligands: Dict[SpecieType, float]  # occupying ligands

    @cached_property
    def total(self):
        return self.free + sum(self.ligands.values())

    def occ(self, spec: SpecieType) -> float:
        return self.ligands[spec]

    def occ_ratio(self, spec: SpecieType) -> float:
        return self.occ(spec) / self.total

    def update(self, ligand: SpecieType, balance: float = 0) -> None:
        shift = self.bind(balance, ligand)
        self.free -= shift
        self.ligands[ligand] += shift

    def bind(self, val: float, st: SpecieType):
        return min(max(val, -self.ligands[st]), self.free)


class Cell:
    specie_receptors: Dict[SpecieType, Occupancy]
    cell_pathway: SimplePathway
    production_rate: float
    timescale: float = 1e-2

    def __init__(self, pathway: SimplePathway, output_rate: float,
                 receptor_count: Dict[SpecieType, int]):
        self.specie_receptors = {
            SpecieType.A: Occupancy(receptor_count[SpecieType.A],
                                    {SpecieType.A: 0, SpecieType.cA: 0}),
            SpecieType.B: Occupancy(receptor_count[SpecieType.B],
                                    {SpecieType.B: 0, SpecieType.cB: 0}),
            SpecieType.C: Occupancy(receptor_count[SpecieType.C],
                                    {SpecieType.C: 0})
        }
        self.specie_receptors[SpecieType.cA] = self.specie_receptors[SpecieType.A]
        self.specie_receptors[SpecieType.cB] = self.specie_receptors[SpecieType.B]
        self.cell_pathway = pathway
        self.production_rate = output_rate

    def update_occupancy(self, comp: Compartment) -> 'Cell':
        for st in self.specie_receptors:
            conc = comp.get_concentration(st)
            occ_balance = (self[st].free * conc * st.kon -
                           self[st].occ(st) * st.koff) * \
                          self.timescale
            self[st].update(st, occ_balance)
            comp.species[st] -= occ_balance
        return self

    def get_occupancy_ratio(self, spec: SpecieType):
        return self[spec].occ_ratio(spec)

    @property
    def occupancy_ratios(self) -> Dict[SpecieType, float]:
        return {st: self.get_occupancy_ratio(st) for st
                in self.specie_receptors.keys()}

    @property
    def production(self) -> float:
        self.cell_pathway.update_inputs(self.occupancy_ratios)
        return self.production_rate * \
            self.cell_pathway \
                .update_inputs(self.occupancy_ratios) \
                .output

    def __getitem__(self, spec: SpecieType) -> Occupancy:
        return self.specie_receptors[spec]

    def __setitem__(self, spec: SpecieType, occ: Occupancy):
        self.specie_receptors[spec] = occ
