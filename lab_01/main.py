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
