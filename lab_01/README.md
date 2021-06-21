# Лабораторная работа №1

**Цель работы:** Получение  навыков  построения  алгоритмаинтерполяции таблично.
заданных функций полиномами Ньютона и Эрмита.

## 1. Исходные данные

1. Таблица функции и её производных

| x    |  y         | y'      |
|------|------------|---------|
| 0.00 |  1.000000  | -1.00000|
| 0.15 |  0.838771  | -1.14944|
| 0.30 |  0.655336  | -1.29552|
| 0.45 |  0.450447  | -1.43497|
| 0.60 |  0.225336  | -1.56464|
| 0.75 | -0.018310  | -1.68164|
| 0.90 | -0.278390  | -1.78333|
| 1.05 | -0.552430  | -1.86742|

2. Степень аппроксимирующего полинома - n.
3. Значение аргумента для которого выполняется интерполяция.

## 2. Код программы

Код программы представлен на листингах 1-2.

### Листинг 1. main.py

```python
import interpolation as interp


def read_table(filename):
    """Считывает таблицу (x, y, производная) из файла"""
    table = []
    try:
        with open(filename, "r") as f:
            for line in f:
                point_data = list(map(float, line.split()))
                point = interp.FuncPoint(point_data[0], point_data[1],
                                         point_data[2])
                table.append(point)
    except FileNotFoundError:
        print("Файл не найден")
        exit()
    except IndexError:
        print("Недостаточно данных в файле")
    except ValueError:
        print("Встречено не число в файле")
        exit()

    return table


def calc_newton(table, polynomial_power, x_arg):
    table_conf = interp.make_configuration(table, polynomial_power, x_arg,
                                           interp.POL_TYPE_NEWTON)
    row_dd = interp.calc_table_divided_differences(table_conf,
                                                   interp.POL_TYPE_NEWTON)[0]
    y_approx = interp.calc_polynomial(table_conf, row_dd, x_arg)
    return y_approx


def calc_hermite(table, polynomial_power, x_arg):
    table_conf = interp.make_configuration(table, polynomial_power, x_arg,
                                           interp.POL_TYPE_HERMITE)
    row_dd = interp.calc_table_divided_differences(table_conf,
                                                   interp.POL_TYPE_HERMITE)[0]
    y_approx = interp.calc_polynomial(table_conf, row_dd, x_arg)
    return y_approx


def main():
    table_filename = "table.txt"
    table = read_table(table_filename)

    x_arg = 0.525

    print(f"x = {x_arg}")
    print("n |     Ньютон |      Эрмит")
    print("--|------------|-----------")
    for n in [0, 1, 2, 3, 4]:
        print("%d | %10.7f | %10.7f" % (n, calc_newton(table, n, x_arg),
                                        calc_hermite(table, n, x_arg)))
    print()

    rev_table = [interp.FuncPoint(point.y, point.x, None) for point in table]
    print("Корень функции при помощи обратной интерполяции")
    print("n |     Ньютон")
    print("--|-----------")
    for n in [0, 1, 2, 3, 4]:
        print("%d | %10.7f" % (n, calc_newton(rev_table, n, 0.0)))


if __name__ == "__main__":
    main()
```

### Листинг 2. interpolation.py

```python
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
```

## 3. Результаты работы

1. Значения y(x) при степенях полиномов Ньютона и Эрмита n=1, 2, 3  и 4
при фиксированном x, например, x=0.525 (середина интервала 0.45-0.60).
Результаты свести в таблицу для сравнения полиномов.

`x = 0.525`

|n |     Ньютон |      Эрмит|
|--|------------|-----------|
|0 |  0.2253360 |  0.2253360|
|1 |  0.3378915 |  0.3426840|
|2 |  0.3402084 |  0.3402877|
|3 |  0.3403138 |  0.3403228|
|4 |  0.3403245 |  0.3403240|

2. Найти корень заданной выше табличной функции с помощью обратной
интерполяции, используя полином Ньютона.

`Корень функции при помощи обратной интерполяции`

|n |     Ньютон|
|--|-----------|
|0 |  0.7500000|
|1 |  0.7387275|
|2 |  0.7390461|
|3 |  0.7390948|
|4 |  0.7390877|

## 4. Вопросы по защите лабораторой работы

1. Будет ли работать программа при степени полинома n=0?\
Программы будет работать, полином будет построен на одном узле.
2. ак практически оценить погрешность интерполяции?\
Почему сложно применить для этих целей теоретическую оценку?\
Погрешность многочлена Ньютона можно оценить по формуле:

![Observational error formula](./readme_imgs/observational_error.svg)

где

![Observational error formula exp 01](./readme_imgs/observational_error_sub_01.svg)

![Observational error formula exp 02](./readme_imgs/observational_error_sub_02.svg)

3. Если в двух точках заданы значения функции и ее первых производных,
то полином какой минимальной степени может быть построен на этих точках?\
Минимум  - 0, максимум - 3.
4. В каком месте алгоритма построения полинома существенна информация
об упорядоченности аргумента функции (возрастает, убывает)?\
Порядок нумерации узлов безразличен.
5. Что такое выравнивающие переменные и как их применить
для повышения точности интерполяции?\
Выравнивающие переменные - это такие переменные η=η(y) ξ=ξ(x), что график
η(ξ) близок к прямой, хотя бы на отдельных участках. В случае
быстроизменяющихся функций интерполяцию проводят на этих переменных, а зачем
приводят обратно к (y, x). Это помогает избегать составления таблиц больших
объемов.
