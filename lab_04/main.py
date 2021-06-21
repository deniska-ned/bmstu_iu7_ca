from sys import argv
from typing import Tuple

import matplotlib.pyplot as plt

from meansquare import *


def read_points(fname: str) -> List[Point]:
    dots = []

    with open(fname) as fin:
        for line in fin.readlines():
            dots += [Point(*list(map(float, line.split()[:3])))]

    return dots


def print_points(dots: List[Point]) -> None:
    print("| {:>8s} | {:>8} | {:>8} |".format("X", "Y", "Weight"))
    print("| {:>8s} | {:>8} | {:>8} |".format('-' * 8, '-' * 8, '-' * 8))
    for i in dots:
        print("| {:8.2f} | {:8.2f} | {:8.2f} |".format(i.x, i.y, i.weight))


def plot(dots: List[Point], approx: List[Tuple[int, List[Point]]]) -> None:
    x, y = [p.x for p in dots], [p.y for p in dots]
    plt.clf()

    plt.title("Meansquare method approximation ")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(which='minor', color='k', linestyle=':')
    plt.grid(which='major', color='k')

    plt.plot(x, y, "o", label = "table points")

    for a in approx:
        plt.plot([p.x for p in a[1]], [p.y for p in a[1]],
                 label="p = {:d}".format(a[0]))

    plt.legend()
    plt.show()


def main():
    dots = read_points(argv[1])

    print("Points from file {:s}".format(argv[1]))
    print_points(dots)

    approxs = []
    for deg in [1, 2, 4, 9]:
        slae = SLAE().build(dots, deg)
        slae = slae.solve()
        approxs.append((deg, Approx().get_coeffs(slae).build(dots)))

    plot(dots, approxs)


if __name__ == "__main__":
    main()
