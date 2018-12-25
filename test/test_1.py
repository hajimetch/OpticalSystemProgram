# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 16:01:37 2017

@author: hajimetch
"""

from typing import Tuple, Dict


class Parts():
    @property
    def radius(self) -> float:
        return self.__radius

    @radius.setter
    def radius(self, input_radius: float) -> None:
        self.__radius = input_radius

    @property
    def r_index(self) -> float:
        return self.__r_index

    @r_index.setter
    def r_index(self, input_r_index: float) -> None:
        self.__r_index = input_r_index

    @property
    def rsrc(self) -> Dict[str, float]:
        return self.__rsrc

    @rsrc.setter
    def rsrc(self, input_rsrc: Dict[str, float]) -> None:
        self.__rsrc = input_rsrc

    def set_rsrc(self, front: float, back: float) -> None:
        self.__rsrc = {'front': front, 'back': back}


if __name__ == '__main__':
    parts = Parts()
    parts.radius = 20
    print(parts.radius)
    parts.set_rsrc(front=2.0, back=3.0)
    print(parts.rsrc)
