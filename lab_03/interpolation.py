import math
from typing import List
from point import Point


POL_TYPE_NEWTON = "Newton"


def make_configuration(nodes: List[Point], pow_pol, x_arg, pol_type):
    """Формируется  конфигурация  из  (pow_pol+1)  узлов,
    по  возможности симметрично расположенных относительно значения  x_arg"""
    table_conf = sorted(nodes,
                        key=lambda point: math.fabs(point.x - x_arg))
    table_conf = table_conf[:pow_pol + 1]
    table_conf.sort(key=lambda point: point.x)
    return table_conf


def calc_dd(config, table_dd, row_i, col_i, pol_type):
    """Вычисляет разделенную разность для таблицы разделенных разностей"""
    dy = (table_dd[row_i][col_i - 1] - table_dd[row_i + 1][col_i - 1])
    dx = config[row_i].x - config[row_i + col_i].x
    return dy / dx


def calc_table_divided_differences(config, pol_type):
    """Формирует таблицу разделенных разностей функции"""
    trc = len(config)
    table_dd = [[0.0] * row_len for row_len in range(trc, 0, -1)]

    for i in range(len(config)):
        table_dd[i][0] = config[i].y

    for col_i in range(1, trc):
        row_count = trc - col_i
        for row_i in range(0, row_count):
            tmp = calc_dd(config, table_dd, row_i, col_i, pol_type)
            table_dd[row_i][col_i] = tmp

    return table_dd


def calc_polynomial(config, row_dd, x_arg):
    """Строит полином Ньютона"""
    res = row_dd[0]
    x_deltas = 1
    for i in range(len(config) - 1):
        x_deltas *= x_arg - config[i].x
        res += x_deltas * row_dd[i + 1]

    return res


def calc_newton(nodes, polynomial_power, x_arg):
    """Вычисляет полином Ньютона"""
    nodes_conf = make_configuration(nodes, polynomial_power, x_arg,
                                    POL_TYPE_NEWTON)
    row_dd = calc_table_divided_differences(nodes_conf, POL_TYPE_NEWTON)[0]
    y_approx = calc_polynomial(nodes_conf, row_dd, x_arg)
    return y_approx


if __name__ == "__main__":
    print("This is the package file")
