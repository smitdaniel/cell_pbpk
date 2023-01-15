from dataclasses import dataclass
from enum import Enum, unique
from functools import cached_property
from typing import Dict


@dataclass
class SpeciesData:
    name: str
    kon: float = 1e-4
    koff: float = 1e-8
    decay: float = 1e-10


@unique
class SpecieType(Enum):
    A = SpeciesData("A")
    B = SpeciesData("B")
    C = SpeciesData("C")
    O = SpeciesData("O")

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

    def __init__(self, volume: float):
        self.volume = volume
        self.species = {st: 0 for st in SpecieType}

    def __getitem__(self, spec: SpecieType) -> float:
        return self.species[spec]

    def __setitem__(self, key, value):
        self.species[key] = value

    def get_concentration(self, spec: SpecieType) -> float:
        return self[spec] / self.volume

    def update_concentration(self) -> None:
        for st in self.species:
            self[st] -= self[st] * st.decay

    def inject(self, a: int = 1000, b: int = 1000, c: int = 1000,
               o: int = 0) -> None:
        for st, update in zip(self.species, [a, b, c, o]):
            self[st] += update
