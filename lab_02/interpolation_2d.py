from typing import List, Tuple, Dict

import interpolation as interp_1d


class TableInterp2d:
    """Хранит таблицу таблично заданной функции 2х переменных"""

    def __init__(self, x_args: List[float], y_args: List[float],
                 values: List[List[float]]):
        assert len(y_args) == len(values)
        assert len(x_args) == len(values[0])

        self.x_args = x_args
        self.y_args = y_args
        self.values = values


def interpolate_2d(table: TableInterp2d, x_pow: int, y_pow: int,
                   x_arg: float, y_arg: float):
    """Выполянет интерполяцию таблично заданных фунций 2х переменных"""

    yy = ({"index": i, "value": y} for i, y in enumerate( table.y_args))
    yy = sorted(yy, key=lambda p: abs(p["value"] - y_arg))
    yy = sorted(yy[:y_pow + 1], key=lambda p: abs(p["index"]))

    nodes_yz: List[interp_1d.FuncPoint] = []

    for y_row_data in yy:
        y_row_value = y_row_data["value"]
        y_row_index = y_row_data["index"]

        nodes_xz: List[interp_1d.FuncPoint] = []

        for i, x in enumerate(table.x_args):
            nodes_xz.append(interp_1d.FuncPoint(x, table.values[y_row_index][i]))

        z_value = interp_1d.calc_newton(nodes_xz, x_pow, x_arg)
        nodes_yz.append(interp_1d.FuncPoint(y_row_value, z_value))

    return interp_1d.calc_newton(nodes_yz, y_pow, y_arg)


if __name__ == "__main__":
    print("This is the package file")
