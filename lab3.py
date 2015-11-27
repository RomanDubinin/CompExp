import copy

N = 4
M = 20 - N
A = 0.1 * M + 0.01 * N
B = 0.2 * M + 0.02 * N + 0.01*M*N + 0.001*N*N

x_star = [
    [1],
    [2],
    [0.1*M + 0.01*N]
]

set_of_equations = [
    [1.2345, 3.1415, 1,         7.5175 + A],
    [2.3456, 5.9690, 0,         14.2836],
    [3.4567, 2.1828, 2 + 0.1*N, 7.8223 + B]
]

ACCURACY_AFTER_DOT = 0


class MethodError(Exception):
    pass


def round(value, accuracy):
    return int((value * (10 ** accuracy))) / (10 ** accuracy)


def get_indexes_of_no_zero_items(list_):
    return set(i for i, element in enumerate(list_)
               if abs(round(element, ACCURACY_AFTER_DOT - 1)) > 0)


def matrix_to_str(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    return '\n'.join(table)


def matrixmult(A, B):
    if len(A[0]) != len(B):
        raise ValueError("Cannot multiply matrices. Incorrect dimensions.")

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(len(B[0]))] for col in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(A[0])):
                C[i][j] += round(A[i][k] * B[k][j], ACCURACY_AFTER_DOT)
    return C


def calculate_cell_by_gauss(i, j, A, B, C):
    sum_ = 0
    for k in range(0, i):
        sum_ += round(B[i][k] * C[k][j], ACCURACY_AFTER_DOT)
    for k in range(i + 1, len(B[1])):
        sum_ += round(B[i][k] * C[k][j], ACCURACY_AFTER_DOT)

    if j > i:  # this element is in C
        C[i][j] = round((A[i][j] - sum_) / B[i][i], ACCURACY_AFTER_DOT)
    else:  # element in B
        B[i][j] = round((A[i][j] - sum_) / C[i][i], ACCURACY_AFTER_DOT)


def calculate_matrix_B_nad_C_by_gauss(A):
    B = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ]
    C = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0]
    ]

    height = len(A)
    weedth = len(A[0])

    try:
        for k in range(0, weedth):
            for i in range(k, height):
                calculate_cell_by_gauss(i, k, A, B, C)
            for j in range(k + 1, weedth):
                calculate_cell_by_gauss(k, j, A, B, C)
    except ZeroDivisionError:
        raise MethodError("Method can not find solution")

    return B, C


def solve_up_triangle_system(matrix, b_column):
    x_column = [[0] for _ in range(len(b_column))]

    for i in range(len(b_column))[::-1]:
        sum = 0
        for k in range(i):
            sum += round(matrix[i][k] * x_column[k][0], ACCURACY_AFTER_DOT)

        for k in range(i+1, len(b_column)):
            sum += round(matrix[i][k] * x_column[k][0], ACCURACY_AFTER_DOT)

        x_column[i][0] = round((b_column[i][0] - sum) /
                               matrix[i][i], ACCURACY_AFTER_DOT)
    return x_column


def solve_any_transformed_system(set_of_equations):
    x_column = []
    for _ in range(len(set_of_equations[0]) - 1):
        x_column.append([0])

    computed_elements_of_x = []
    num_of_computed_x_on_last_iteartion = 0
    no_zero_indexes_on_last_iteration = []

    while len(computed_elements_of_x) != len(x_column):
        for equation in set_of_equations:
            no_zero_elements_in_current_row = get_indexes_of_no_zero_items(equation[:-1])
            new_no_zero_elements = list(
                no_zero_elements_in_current_row -
                set(no_zero_indexes_on_last_iteration)
            )
            if len(new_no_zero_elements) != 1:
                continue

            j = new_no_zero_elements[0]

            sum_ = 0
            for k in range(0, j):
                sum_ += round(equation[k] * x_column[k][0], ACCURACY_AFTER_DOT)
            for k in range(j + 1, len(set_of_equations[0]) - 1):
                sum_ += round(equation[k] * x_column[k][0], ACCURACY_AFTER_DOT)

            x_column[j][0] = round((equation[-1] - sum_) / equation[j],
                                   ACCURACY_AFTER_DOT)

            computed_elements_of_x.append(j)
            num_of_computed_x_on_last_iteartion += 1
            no_zero_indexes_on_last_iteration = no_zero_elements_in_current_row
            break

        if (num_of_computed_x_on_last_iteartion ==
           len(computed_elements_of_x) - 1):
            raise MethodError("Method is unable to find solution!")
    return x_column


def compact_gauss_scheme(set_of_equations):
    B, C = calculate_matrix_B_nad_C_by_gauss(set_of_equations)

    y_column = [[row[-1]] for row in C]

    C = [line[:-1] for line in C]

    x_by_method = solve_up_triangle_system(C, y_column)

    print("Matrix B:\n" + matrix_to_str(B) + "\n")
    print("Matrix C:\n" + matrix_to_str(C) + "\n")
    print("column y:\n" + matrix_to_str(y_column) + "\n")
    return x_by_method


def find_max_cell(matrix, start_string_num):
    max_val = matrix[start_string_num][0]
    res_i = start_string_num
    res_j = 0
    for i in range(start_string_num, len(matrix)):
        for j in range(0, len(matrix[i])):
            if abs(matrix[i][j]) > max_val:
                max_val = matrix[i][j]
                res_i = i
                res_j = j
    return res_i, res_j


def swap_matrix_lines(matrix, string_num_1, string_num_2):
    temp = matrix[string_num_1]
    matrix[string_num_1] = matrix[string_num_2]
    matrix[string_num_2] = temp


def do_forward_step_of_gauss_by_main_element(matrix):
    set_of_equations = copy.deepcopy(matrix)

    for start_string_index in range(len(set_of_equations)):
        A = [line[:-1] for line in set_of_equations]

        main_i, main_j = find_max_cell(A, start_string_index)
        swap_matrix_lines(set_of_equations, start_string_index, main_i)

        for equation in set_of_equations[start_string_index + 1:]:
            multiplicator = round(-equation[main_j] /
                                  set_of_equations[start_string_index][main_j],
                                  ACCURACY_AFTER_DOT)
            for j in range(len(equation)):
                equation[j] = round(
                    equation[j] +
                    set_of_equations[start_string_index][j]*multiplicator,
                    ACCURACY_AFTER_DOT
                )
    return set_of_equations


def main_item_gauss(matrix):
    set_of_equations = do_forward_step_of_gauss_by_main_element(matrix)
    print("Set of equations after forward step of method:\n" +
          matrix_to_str(set_of_equations) + "\n")
    return solve_any_transformed_system(set_of_equations)


def do_method_investigation(method, accuracy_list):
    for accuracy in accuracy_list:
        global ACCURACY_AFTER_DOT
        ACCURACY_AFTER_DOT = accuracy

        copy_of_system = []

        for i in range(len(set_of_equations)):
            copy_of_system.append([])
            for j in range(len(set_of_equations[i])):
                copy_of_system[i].append(
                    round(set_of_equations[i][j], ACCURACY_AFTER_DOT))

        try:
            print("\nAccuracy: ", accuracy)
            x_by_method = method(copy_of_system)
            print("Answer:\n" + matrix_to_str(x_by_method) + "\n")
            
        except MethodError as e:
            print(e)


def main():
    import subprocess
    subprocess.call("clear", shell=True)
    subprocess.call("cls", shell=True)

    print("[I]\t\t\tReal answer(X*):\n\n" + matrix_to_str(x_star) + "\n")
    print("=" * 240)

    print("[II]\t\t\tCompact Gauss Scheme: ")
    do_method_investigation(compact_gauss_scheme, accuracy_list=[2, 4, 6])
    print("=" * 240)

    print("[III]\t\t\tGauss Scheme with main element: ")
    do_method_investigation(main_item_gauss, accuracy_list=[2, 4, 6])


if __name__ == '__main__':
    main()
