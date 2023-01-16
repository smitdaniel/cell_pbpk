from dataclasses import dataclass
from enum import Enum, unique
from functools import cached_property
from typing import Dict, Optional


@dataclass
class SpeciesData:
    name: str
    kon: float
    koff: float
    decay: float


@unique
class SpecieType(Enum):
    A = SpeciesData("A", kon=1e-4, koff=1e-4, decay=1e-2)
    cA = SpeciesData("cA", kon=1e-5, koff=1e-4, decay=1e-3)
    B = SpeciesData("B", kon=1e-3, koff=1e-1, decay=1e-3)
    cB = SpeciesData("cB", kon=1e-5, koff=1e-5, decay=1e-6)
    C = SpeciesData("C", kon=1e-3, koff=1e-8, decay=1e-15)

    def __str__(self):
        return self.value.name

    @cached_property
    def kon(self):
        return self.value.kon

    @cached_property
    def koff(self):
        return self.value.koff

    @cached_property
    def decay(self):
        return self.value.decay


class Compartment:
    volume: float
    species: Dict[SpecieType, int]
    injection: Dict[SpecieType, int]

    def __init__(self, injection: Dict[SpecieType, int], volume: float = 1.0):
        self.volume = volume
        self.species = {st: 0 for st in SpecieType}
        self.injection = injection

    def __getitem__(self, spec: SpecieType) -> float:
        return self.species[spec]

    def __setitem__(self, key, value):
        self.species[key] = value

    def get_concentration(self, spec: SpecieType) -> float:
        return self[spec] / self.volume

    def update_species_decay(self) -> None:
        for st in self.species:
            self[st] -= self[st] * st.decay

    def inject(self, inj_dict: Optional[Dict[SpecieType, int]] = None) -> None:
        inj_dict = inj_dict if inj_dict is not None else self.injection
        for st in self.species:
            self[st] += inj_dict.get(st, 0)
