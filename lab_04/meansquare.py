from __future__ import annotations

import numpy as np
from typing import List


class Point:
    def __init__(self, _x: float, _y: float, _w: float):
        self.x = _x
        self.y = _y
        self.weight = _w


class SLAE:
    mat: List[List[float]]
    n: int

    def build(self, ds: List[Point], _n: int) -> SLAE:
        self.n = _n
        self.mat = [[0 for _ in range(self.n + 2)] for _ in range(self.n + 1)]

        for i in range(self.n + 1):
            for j in range(self.n + 1):
                slae_coeffs = 0.0
                expanded_coeff = 0.0
                for k in range(len(ds)):
                    slae_coeffs += (ds[k].weight
                                    * (ds[k].x ** i)
                                    * (ds[k].x ** j))
                    expanded_coeff += ds[k].weight * ds[k].y * (ds[k].x ** i)

                self.mat[i][j] = slae_coeffs
                self.mat[i][self.n + 1] = expanded_coeff

        return self

    def solve(self) -> List[List[float]]:
        for i in range(self.n + 1):
            for j in range(self.n + 1):
                if i == j:
                    continue

                sub_coeff = self.mat[j][i] / self.mat[i][i]
                for k in range(self.n + 2):
                    self.mat[j][k] -= sub_coeff * self.mat[i][k]

        for i in range(self.n + 1):
            divider = self.mat[i][i]
            for j in range(self.n + 2):
                self.mat[i][j] /= divider

        return self.mat


class Approx:
    def __init__(self):
        self.coeffs = []

    def get_coeffs(self, mat: List[List[float]]) -> Approx:
        self.coeffs = [mat[i][len(mat)] for i in range(len(mat))]

        return self

    def build(self, ds: List[Point]) -> List[Point]:
        dots = []

        for i in np.arange(ds[0].x, ds[-1].x + 0.1, 0.1):
            d = Point(i, 0, 0)

            for j in range(len(self.coeffs)):
                d.y += d.x ** j * self.coeffs[j]

            dots += [d]

        return dots
