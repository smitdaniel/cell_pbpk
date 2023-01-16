from typing import Dict, Tuple, Union

from src.compartment import SpecieType


class SimplePathway:
    inputs: Dict[SpecieType, bool]
    threshold: Dict[SpecieType, float]

    def __init__(self, thresh: Dict[SpecieType, float]):
        self.threshold = thresh
        self.inputs = {
            SpecieType.A: False,
            SpecieType.B: False,
            SpecieType.C: False}

    def __getitem__(self,
                    lbl: Union[SpecieType, Tuple[SpecieType, SpecieType]]
                    ) -> bool:
        if isinstance(lbl, SpecieType):
            return self.inputs[lbl]
        elif isinstance(lbl, tuple):
            if lbl == (SpecieType.A, SpecieType.B):
                return self.node_ab
            if lbl == (SpecieType.B, SpecieType.C):
                return self.node_bc

    @property
    def node_ab(self) -> bool:
        return self[SpecieType.A] and not self[SpecieType.B]

    @property
    def node_bc(self) -> bool:
        return self[SpecieType.B] or not self[SpecieType.C]

    def update_inputs(self,
                      occ_ratios: Dict[SpecieType, float]
                      ) -> 'SimplePathway':
        for st in self.inputs:
            self.inputs[st] = self.filter(st, occ_ratios[st])
        return self

    def filter(self, spec: SpecieType, occ_ratio: float) -> bool:
        return occ_ratio >= self.threshold[spec]

    @property
    def output(self) -> bool:
        return self.node_ab or not self.node_bc
