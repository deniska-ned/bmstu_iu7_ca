# Лабораторная работа №2

**Teма:** Построение и программная реализация алгоритма многомерной
интерполяции табличных функций.

**Цель работы:** Получение навыков построения алгоритма интерполяции таблично
заданных фунций двух переменных.

## 1. Исходные данные

1. Таблица функции с количеством узлов 5x5.

| y\x |  0 |  1 |  2 |  3 |  4 |
|-----|----|----|----|----|----|
|   0 |  0 |  1 |  4 |  9 | 16 |
|   1 |  1 |  2 |  5 | 10 | 17 |
|   2 |  4 |  5 |  8 | 13 | 20 |
|   3 |  9 | 10 | 13 | 18 | 25 |
|   4 | 16 | 17 | 20 | 25 | 32 |

2. Степень аппроксимирующих полиномов - n\_x, n\_y.
3. Значение аргументов x, y для которого выполняется интерполяция.

## 2. Код программы

Код программы представлен на листингах 1-3.

### Листинг 1. interpolation.py

```python
import math
from typing import List


POL_TYPE_NEWTON = "Newton"
POL_TYPE_HERMITE = "Hermite"


class FuncPoint:

    def __init__(self, x=0.0, y=0.0, d=None):
        self.x = x
        self.y = y
        self.d = d

    def __repr__(self):
        return f"x = {self.x}; y = {self.y}; d = {self.d}"


def make_configuration(nodes: List[FuncPoint], pow_pol, x_arg, pol_type):
    """Формируется  конфигурация  из  (pow_pol+1)  узлов,
    по  возможности симметрично расположенных относительно значения  x_arg"""
    if pol_type == POL_TYPE_HERMITE:
        nodes = nodes * 2

    table_conf = sorted(nodes,
                        key=lambda point: math.fabs(point.x - x_arg))
    table_conf = table_conf[:pow_pol + 1]
    table_conf.sort(key=lambda point: point.x)
    return table_conf


def calc_dd(config, table_dd, row_i, col_i, pol_type):
    """Вычисляет разделенную разность для таблицы разделенных разностей"""
    if (pol_type == POL_TYPE_HERMITE and col_i == 1 and
            math.fabs(config[row_i].x - config[row_i + 1].x) < 1e-5):
        return config[row_i].d

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


def calc_hermite(nodes, polynomial_power, x_arg):
    """Вычисляет полином Эрмита"""
    nodes_conf = make_configuration(nodes, polynomial_power, x_arg,
                                    POL_TYPE_HERMITE)
    row_dd = calc_table_divided_differences(nodes_conf, POL_TYPE_HERMITE)[0]
    y_approx = calc_polynomial(nodes_conf, row_dd, x_arg)
    return y_approx


if __name__ == "__main__":
    print("This is the package file")
```

### Листинг 2. interpolation_2d.py

```python
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
```

### Листинг 3. main.py

```python
import interpolation_2d as interp_2d


def main():
    x_args = [0, 1, 2, 3, 4]
    y_args = [0, 1, 2, 3, 4]
    values = [[0, 1, 4, 9, 16],
              [1, 2, 5, 10, 17],
              [4, 5, 8, 13, 20],
              [9, 10, 13, 18, 25],
              [16, 17, 20, 25, 32]
              ]

    table = interp_2d.TableInterp2d(x_args, y_args, values)

    x_arg = 1.5
    y_arg = 1.5

    print(f"x = {x_arg}")
    print(f"y = {y_arg}")
    print()

    print("| %9s | x_pow = %d | x_pow = %d | x_pow = %d |" % ("", 1, 2, 3))
    print("|%11s|%11s|%11s|%11s|" % ("-" * 11, "-" * 11, "-" * 11, "-" * 11))

    for y_pow in [1, 2, 3]:
        print("| y_pow = %d" % y_pow, end="")
        for x_pow in [1, 2, 3]:
            z_value = interp_2d.interpolate_2d(table, x_pow, y_pow,
                                               x_arg, y_arg)
            print(" | %9.2lf" % z_value, end="")
        print(" |")


if __name__ == "__main__":
    main()
```

## 3. Результаты работы

Результат интерполяции z(x, y) при степенях полиномов 1, 2, 3 для `x = 1.5`,
`y = 1.5`.

|            | x\_pow = 1 | x\_pow = 2 | x\_pow = 3 |
|------------|------------|------------|------------|
| y\_pow = 1 |       5.00 |       4.75 |       4.75 |
| y\_pow = 2 |       4.75 |       4.50 |       4.50 |
| y\_pow = 3 |       4.75 |       4.50 |       4.50 |

## 4. Вопросы по защите лабораторой работы

1. Пусть производящая функция таблицы суть `z(x, y)=x^2 + y^2`. Область
определения по x и y `0-5` и `0-5`. Шаги по переменных равны 1. Степени
n\_x = n\_y = 1, x = y = 1.5. Приведите по шагам те значения фунции, которые
получаются в ходе последовательных интерполяций по строкам и столбцу.\
По строкам:

|1.0|3.5|
|-|-|
|2.0|5.0|

По столбцу: 5.0

2. Какова миниамльная степень двумерного полинома, построенного на четырех
узлах? На шести узлах?\
На четырех: 3; на шести: 5. Так как количество узлов для полинома n степени:
n + 1.

3. Предложите алгоритм двумерно интерполяции при хаотичном расположении узлов,
т.е. когда таблицы функции на регульрной сетке нет, и метод последовательной
интерполяции не работает. Какие имеются ограничения на расположение узлов
при разных степенях полинома?\
Ограничиваясь интерполяционным полиномом первой степени имеем
`z = a + b * x + c * y`, и его коээфициенты находят по 3м ближайшим узлам
`z_i = a + b * x_i + c * y_i`.\
Ограчение в расположении: при полиноме первой степени узлы не могут лежать на
одной прямой, при 2ой степени: на одной плоскости и т.д.

4. Пусть на каком-либо языке программирования написана функция, выполняющая
интерполяцию по двум переменных. Опишите алгоритм использования этой
функции для интерполяции по трем переменных.\
    1. Выбрать отрезки для интерполяции по 3м измерениям (по x, y, z).
    2. Использование двумерной интерполяции в каждой из плоскостей
       при фиксирвоанных z.
    3. Произвести одномерную интерполяцию по полученным данным.

5. Можно ли при последовательной интерполяции по разным направления
использовать полиномы несовпадающих стпеней или даже разные методы одномерной
интерполяции, например, полином Ньютона и сплайн?\
Можно.

6. Опишите алгоритм двумерной интерполяции на тругольной конфигурации узлов.\
При треугольной конфигурации расположения узлов степень многочлена будет
минимальной. Многочлен n-й степени в форме Ньютона для двумерной интерполяции
в этом случае можно представить как обобщение одномерного варианта записи:
<!-- P_n(x, y) = \sum_{i=0}^{n} \sum_{j=0}^{n-1} z(x_0, \cdots , x_i, y_0, \cdots , y_i) \prod_{p=0}^{i-1} (x - x_p) \prod_{q=0}^{j-1} (y - y_p) -->

![two-dimensional interpolation polynomial](./readme_imgs/two-dimensional_interpolation_polynomial.svg)
