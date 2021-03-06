# Лабораторная работа №6

**Тема:** Построение и программная реализация алгоритмов численного
дифференцирования.

**Цель работы:** Получение навыков построения алгоритма вычисления производных
от сеточных функций.

## 1. Исходные данные

Задана табличная (сеточная) функция. Имеется информация, что закономерность,
представленная этой таблицей, может быть описана формулой

![](./images/equation_01.svg)

| x | y     | 1 | 2 | 3 | 4 | 5 |
|---|-------|---|---|---|---|---|
| 1 | 0.571 |   |   |   |   |   |
| 2 | 0.889 |   |   |   |   |   |
| 3 | 1.091 |   |   |   |   |   |
| 4 | 1.231 |   |   |   |   |   |
| 5 | 1.333 |   |   |   |   |   |
| 6 | 1.412 |   |   |   |   |   |

Вычислить первые разностные производные от функции и занести их в столбцы
(1)-(4) таблицы:

1. односторонняя разностная производная ,
2. центральная разностная производная,
3. 2-я формула Рунге с использованием односторонней производной,
4. введены выравнивающие переменные.

В столбец 5 занести вторую разностную производную.

## 2. Код программы

### Листинг 1. main.py

```python
from functions import *


def main():
    table = [[1.0, 0.571],
             [2.0, 0.889],
             [3.0, 1.091],
             [4.0, 1.231],
             [5.0, 1.333],
             [6.0, 1.412]]

    left_diff(table)
    centre_diff(table)
    second_Runge(table)
    align_vars(table)
    second_der(table)

    print('|   x   |   y   |  left | centre| Runge | align | second|')
    print('|-------|-------|-------|-------|-------|-------|-------|')
    for string in table:
        print('|', end='')
        for field in string:
            if field == '   -   ':
                print(field + '|', end='')
            else:
                print('{:7.4f}|'.format(field), end='')
        print()
    

if __name__ == "__main__":
    main()
```

### Листинг 2. functions.py

```python
def left_diff_formula(dot1, dot2):
    return ((dot1[1] - dot2[1]) /
            (dot1[0] - dot2[0]))


def left_diff(table):
    table[0].append('   -   ')
    for i in range(1, len(table)):
        table[i].append(left_diff_formula(table[i], table[i - 1]))

    
def centre_diff(table):
    table[0].append('   -   ')
    for i in range(1, len(table) - 1):
        table[i].append((table[i + 1][1] - table[i - 1][1]) /
                        (table[i + 1][0] - table[i - 1][0]))
    table[i + 1].append('   -   ')


def second_Runge(table):
    table[0].append('   -   ')
    table[1].append('   -   ')
    for i in range(2, len(table)):
        table[i].append(left_diff_formula(table[i], table[i - 1]) * 2 -
                        left_diff_formula(table[i], table[i - 2]))


def align_vars(table):
    new_table = []
    for dot in table:
        new_table.append([1 / dot[0], 1 / dot[1]])

    table[0].append('   -   ')
    for i in range(1, len(new_table)):
        table[i].append(left_diff_formula(new_table[i], new_table[i - 1]) *
                        table[i][1] ** 2 / table[i][0] ** 2)


def second_der(table):
    table[0].append('   -   ')
    for i in range(1, len(table) - 1):
        table[i].append((table[i + 1][1] + table[i - 1][1] - table[i][1] * 2) /
                        (table[i + 1][0] - table[i][0]) ** 2)
    table[i + 1].append('   -   ')
```

## 3. Результаты работы

Заполненная таблица с краткими комментариями по поводу использованных
формул и их точности

|   x   |   y   |  left | centre| Runge | align | second|
|-------|-------|-------|-------|-------|-------|-------|
| 1.0000| 0.5710|   -   |   -   |   -   |   -   |   -   |
| 2.0000| 0.8890| 0.3180| 0.2600|   -   | 0.2475|-0.1160|
| 3.0000| 1.0910| 0.2020| 0.1710| 0.1440| 0.1653|-0.0620|
| 4.0000| 1.2310| 0.1400| 0.1210| 0.1090| 0.1185|-0.0380|
| 5.0000| 1.3330| 0.1020| 0.0905| 0.0830| 0.0884|-0.0230|
| 6.0000| 1.4120| 0.0790|   -   | 0.0675| 0.0697|   -   |

### 3.1 Левая разностная производная

![](./images/equation_02.svg)

Получается из разложения функции в ряд Тейлора

Точность: первый порядок точности относительно шага h.

### 3.2 Центральная разностная производная

![](./images/equation_03.svg)

Получается вычитаем разложения ф-и в ряд Тейлора для `y_n+1` из
разложения ф-и в ряд Тейлора для `y_n-1`

Точность: второй порядок точности относительно шага h.

### 3.3 2-я формула Рунге с использованием односторонней производной:

формула Рунге для левой разностной производной

![](./images/equation_04.svg)

где

![](./images/equation_05.svg)

Точность: 2

### 3.4 Введение выравнивающих переменных:

Выравниваниющие переменные вводятся для преобразования исходной
линии в прямую. Введем следующие выравнивающие переменные:

![](./images/equation_06.svg)

![](./images/equation_07.svg)

Зависимость принимает вид:

![](./images/equation_08.svg)

Формула приобретает вид

![](./images/equation_09.svg)

где

![](./images/equation_10.svg)

Точность: абсолютная

### 3.5 Вторая разностная производная:

![](./images/equation_11.svg)

Получается вычитаем разложения ф-и в ряд Тейлора для `y_n+1` из
разложения ф-и в ряд Тейлора для `y_n-1`

Точность: второй порядок точности относительно шага h.

## 4. Вопросы по защите лабораторой работы

1. Получить формулу порядка точности `O(h^2)` для первой разностной
производной `y'_N` в крайнем правом узле `x_N`.

![](./images/equation_12.svg)

Исключаем слагаемое содержащее `h^2`, тем самым получим трехчленную
формулу

![](./images/equation_13.svg)

2. Получить формулу порядка точности `O(h^2)` для второй разностной производной
`y''_0` в крайнем левом узле `x_0`.

![](./images/equation_14.svg)

![](./images/equation_15.svg)

![](./images/equation_16.svg)

3. Используя 2-ую формулу Рунге, дать вывод выражения (9) из Лекции 7 для первой
производной `y'_0` в левом крайнем узле

![](./images/equation_17.svg)

Формула Рунге:

![](./images/equation_18.svg)

![](./images/equation_19.svg)

4. Любым способом из Лекций 7, 8 получить формулу порядка точности `O (h^3)`
для первой разностной производной `y'_0` в крайнем левом узле `x_0`.

![](./images/equation_20.svg)

Исключаем слагаемое содержащее `h^2`, тем самым получим трехчленную
формулу

![](./images/equation_21.svg)
