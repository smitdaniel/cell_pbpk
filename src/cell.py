from dataclasses import dataclass
from functools import cached_property
from typing import Dict, Union

from src.compartment import SpecieType, Compartment
from src.pathway import Pathway, SimplePathway


@dataclass
class Occupancy:
    free: int
    occ: int = 0

    @cached_property
    def total(self):
        return self.free + self.occ

    @property
    def ratio(self) -> float:
        return self.occ / self.total

    def update_occ_balance(self, balance: float, occ_name: str = 'occ') -> None:
        balance = round(balance)
        self.free = self.bind(self.free - balance)
        setattr(self, occ_name, self.bind(getattr(self, occ_name) + balance))

    def bind(self, val: int):
        return min(max(val, 0), self.total)


@dataclass
class CompetitiveOccupancy(Occupancy):
    cmp_occ: int = 0

    @cached_property
    def total(self):
        return super().total + self.cmp_occ


class Cell:
    specie_receptors: Dict[SpecieType, Union[Occupancy, CompetitiveOccupancy]]
    cell_pathway: Pathway
    production_rate: float

    def __init__(self, pathway: Pathway, output_rate: float):
        self.specie_receptors = {
            SpecieType.A: Occupancy(100),
            SpecieType.B: CompetitiveOccupancy(500),
            SpecieType.C: Occupancy(5000)
        }
        self.cell_pathway = pathway
        self.production_rate = output_rate

    def update_occupancy(self, comp: Compartment) -> 'Cell':
        for st in self.specie_receptors:
            conc = comp.get_concentration(st)
            occ_balance = self[st].free * conc * st.kon - self[st].occ * st.koff
            self[st].update_occ_balance(occ_balance)
        conc = comp.get_concentration(SpecieType.O)
        cmp_balance = self[SpecieType.B].free * conc * SpecieType.O.kon - \
                      self[SpecieType.B].cmp_occ * SpecieType.O.koff
        self[SpecieType.B].update_occ_balance(cmp_balance, 'cmp_occ')
        return self

    def get_occupancy_ratio(self, spec: SpecieType):
        return self[spec].ratio

    @property
    def occupancy_ratios(self) -> Dict[SpecieType, float]:
        return {st: self.get_occupancy_ratio(st) for st
                in self.specie_receptors.keys()}

    @property
    def production(self) -> int:
        return round(self.production_rate *
                     self.cell_pathway \
                     .update_inputs(self.occupancy_ratios) \
                     .output
                     )

    def __getitem__(self, spec: SpecieType) -> Occupancy:
        return self.specie_receptors[spec]

    def __setitem__(self, spec: SpecieType, occ: Occupancy):
        self.specie_receptors[spec] = occ


if __name__ == "__main__":
    cmp = Compartment(1)
    cmp.inject()

    pathway = SimplePathway({
        SpecieType.A: 0.5,
        SpecieType.B: 0.7,
        SpecieType.C: 0.3})
    c1 = Cell(pathway, 10)
    produced = c1.update_occupancy(cmp).production
    cmp.inject(0, 0, 0, produced)
    pass
