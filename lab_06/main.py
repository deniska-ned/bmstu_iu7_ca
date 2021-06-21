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
