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
