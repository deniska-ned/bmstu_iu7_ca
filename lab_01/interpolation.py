import copy
import math


POL_TYPE_NEWTON = "Newton"
POL_TYPE_HERMITE = "Hermite"


class FuncPoint:

    def __init__(self, x=0.0, y=0.0, d=None):
        self.x = x
        self.y = y
        self.d = d

    def __repr__(self):
        return f"x = {self.x}; y = {self.y}; d = {self.d}"


def make_configuration(table, pow_pol, x_arg, pol_type):
    """Формируется  конфигурация  из  (pow_pol+1)  узлов,
    по  возможности симметрично расположенных относительно значения  x_arg"""
    if pol_type == POL_TYPE_HERMITE:
        table = table * 2

    table_conf = sorted(table,
                        key=lambda point: math.fabs(point.x - x_arg))
    table_conf = table_conf[:pow_pol + 1]
    table_conf.sort(key=lambda point: point.x)
    return table_conf


def calc_dd(table_conf, table_dd, row_i, col_i, pol_type):
    """Вычисляет раделенную разность для таблицы разделенных разностей"""
    if (pol_type == POL_TYPE_HERMITE and col_i == 1 and
            math.fabs(table_conf[row_i].x - table_conf[row_i + 1].x) < 1e-5):
        return table_conf[row_i].d

    dy = (table_dd[row_i][col_i - 1] - table_dd[row_i + 1][col_i - 1])
    dx = table_conf[row_i].x - table_conf[row_i + col_i].x
    return dy / dx


def calc_table_divided_differences(table_conf, pol_type):
    """Формирует таблицу разделенных разностей функции"""
    trc = len(table_conf)
    table_dd = [[0.0] * row_len for row_len in range(trc, 0, -1)]

    for i in range(len(table_conf)):
        table_dd[i][0] = table_conf[i].y

    for col_i in range(1, trc):
        row_count = trc - col_i
        for row_i in range(0, row_count):
            tmp = calc_dd(table_conf, table_dd, row_i, col_i, pol_type)
            table_dd[row_i][col_i] = tmp

    return table_dd


def calc_polynomial(table_conf, row_dd, x_arg):
    """Строит полином Ньютона"""
    res = row_dd[0]
    x_deltas = 1
    for i in range(len(table_conf) - 1):
        x_deltas *= x_arg - table_conf[i].x
        res += x_deltas * row_dd[i + 1]

    return res


if __name__ == "__main__":
    print("This is the package file")
