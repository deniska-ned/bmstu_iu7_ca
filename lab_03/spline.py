from __future__ import annotations
from typing import List

from point import Point


class Spline:
    def __init__(self, _points: List[Point]) -> None:
        self.points = _points

    def get_pos(self, d: Point) -> int:
        i = 1

        while i < len(self.points) and self.points[i].x < d.x:
            i += 1

        return i - 1

    def calc(self, x: float) -> Point:
        x_args = [d.x for d in self.points]
        y_args = [d.y for d in self.points]

        a = y_args[:-1]

        c = [0] * (len(self.points) - 1)
        ksi_coef, eta_coef = [0, 0], [0, 0]

        for i in range(2, len(self.points)):
            xhi = x_args[i] - x_args[i - 1]
            xhi_1 = x_args[i - 1] - x_args[i - 2]
            yhi = y_args[i] - y_args[i - 1]
            yhi_1 = y_args[i - 1] - y_args[i - 2]

            fi = 3 * (yhi / xhi - yhi_1 / xhi_1)

            ksi_coef.append(-xhi
                            / (xhi_1 * ksi_coef[i - 1] + 2 * (xhi_1 + xhi)))
            eta_coef.append((fi - xhi_1 * eta_coef[i - 1])
                            / (xhi_1 * ksi_coef[i - 1] + 2 * (xhi_1 + xhi)))

        c[len(self.points) - 2] = eta_coef[-1]

        for i in range(len(self.points) - 2, 0, -1):
            c[i - 1] = ksi_coef[i] * c[i] + eta_coef[i]

        b, d = [], []
        for i in range(1, len(self.points) - 1):
            xhi = x_args[i] - x_args[i - 1]
            yhi = y_args[i] - y_args[i - 1]
            b.append(yhi / xhi - (xhi * (c[i] + 2 * c[i - 1])) / 3)
            d.append((c[i] - c[i - 1]) / (3 * xhi))

        b.append((y_args[-1] - y_args[-2]) / (x_args[-1] - x_args[-2]) -
                 ((x_args[-1] - x_args[-2]) * 2 * c[-1]) / 3)
        d.append(-c[len(self.points) - 2] / (3 * (x_args[-1] - x_args[-2])))

        pos = self.get_pos(Point(x, 0))

        res = (a[pos]
               + b[pos] * (x - self.points[pos].x)
               + c[pos] * (x - self.points[pos].x) ** 2
               + d[pos] * (x - self.points[pos].x) ** 3)

        return Point(x, res)


if __name__ == "__main__":
    print("This is package file")
