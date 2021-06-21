from typing import List
from spline import Point

from spline import Spline

import interpolation


def print_points(points: List[Point]) -> None:
    print("| {:>8} | {:>8} |".format("x", "y"))
    print("|----------|----------|")
    for i in points:
        print("| {:8.2f} | {:8.2f} |".format(i.x, i.y))


def print_table(x_values, y_values):
    print("|        x | spline y | newton y |")
    print("|----------|----------|----------|")

    assert len(x_values) == len(y_values)
    for i in range(len(x_values)):
        print("| {:>8.4f} | {:>8.4f} | {:>8.4f} |".format(
            x_values[i], y_values[i]['spline'].y, y_values[i]['newton'].y))


def main():
    points = [
        Point(0, 0),
        Point(1, 1),
        Point(2, 4),
        Point(3, 9),
        Point(4, 16),
        Point(5, 25),
        Point(6, 36),
        Point(7, 49),
        Point(8, 64),
        Point(9, 81),
        Point(10, 100)
    ]

    print_points(points)
    print()

    x_values = [0.5, 5.5]
    y_values = [{'spline': None, 'newton': None}] * len(x_values)

    for i, x in enumerate(x_values):
        spline_point = Spline(points).calc(x)

        newton_point = Point(x, interpolation.calc_newton(points, 3, x))

        y_values[i] = {'spline': spline_point, 'newton': newton_point}

    print_table(x_values, y_values)


if __name__ == "__main__":
    main()
