# Лабораторная работа №3

**Teма:** Построение и программная реализация алгоритма сплайн-интерполяции
табличных функций.

**Цель работы:** Получение навыков владения методоами интерполяции таблично 
заданных функций с помощью кубических сплайнов.

## 1. Исходные данные

`y = x^2`

|  x |  y  |
|----|-----|
|  0 |   0 |
|  1 |   1 |
|  2 |   4 |
|  3 |   9 |
|  4 |  16 |
|  5 |  25 |
|  6 |  36 |
|  7 |  49 |
|  8 |  64 |
|  9 |  81 |
| 10 | 100 |

Значения x для сравнения: 0.5, 5.5

## 2. Код программы

### Листинг 1. main.py

```python
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
```

### Листинг 2. point.py

```python
class Point:

    def __init__(self, _x: float, _y: float) -> None:
        self.x = _x
        self.y = _y

    def __lt__(self, other):
        return self.x < other.x

    def __repr__(self):
        return "(x = {:.4f}; y = {:.4f})".format(self.x, self.y)


if __name__ == "__main__":
    print("This is package file")
```

### Листинг 3. imterpolation.py

```python
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
```

### Листинг 4. spline.py

```python
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
```


## 3. Результаты работы

|        x | spline y | newton y |
|----------|----------|----------|
|   0.5000 |   0.3415 |   0.2500 |
|   5.5000 |  30.2503 |  30.2500 |

## 4. Вопросы по защите лабораторой работы

1. Получить выражения для коэффициентов кубического сплайна, построенного
на двух точках.

```math
a = x_0
b = (y1 - y0) / (x1 - x0)
c = 0
d = 0
```

2. Выписать все условия для определения коэффициентов сплайна, построенного
на 3-х точках

* Значения 1-го многочлена и интерполируемой функции в 0 и 1 точках совпадают
* Значения 2-го многочлена и интерполируемой функции в 1 и 2 точках совпадают
* Равенство  в  точке  1  первой  и  второй  производных,
  вычисляемых  по коэффициентам на соседних участках
* Второй производной в точках 0 и 2 равно нулю

3. Определить начальные значения прогоночных коэффициентов, если принять,
что для коэффициентов сплайна справедливо C\_1 = C\_2.

```math
c_{i−1} = ξ_i * c_i + η_i => c_1 = ξ_2 * c_2 + η_2 => ξ_2 = 1; η = 0 
```

4. Написать формулу для определения последнего коэффициента сплайна С\_N, чтобы
можно было выполнить обратный ход метода прогонки, если в качестве граничного
условия задано `k * C_{N-1} +m * C_N = p`, где k, m и p - заданные числа.

```math
c_{N−1} = ξ_N * c_N + η_N

c_N = (p - k * η_N) / (k * ξ_N + m)
```
