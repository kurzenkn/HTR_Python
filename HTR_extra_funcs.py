import numpy as np
import pandas as pd

def get_full_matrix(R):
    """Define matrix A for the linear system of equations."""

    #   Create matrices for pressure and flow rate contributions to A
    M_P = np.zeros([31,31])
    M_Q = np.zeros([31,31])

    #   Define the coefficient contributions to matrix A from known pressure values
    M_P[0,0] = 1
    M_P[1,4] = 1
    M_P[2,12] = 1

    #   Define the coefficient contributions to matrix A from the Hagen-Poiseulle Relations
    p_plus_i = [0, 1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 8, 9, 10, 11]
    p_neg_i = [1, 2, 3, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12, 12]
    R_index = range(0, 18)
    HP_eq_num = range(3, 21)
    for i in range(0,len(HP_eq_num)):
        M_P[HP_eq_num[i], p_plus_i[i]] = 1
        M_P[HP_eq_num[i], p_neg_i[i]] = -1
        M_Q[HP_eq_num[i], R_index[i] + 13] = -R[0, R_index[i]]

    #   Define the coefficient contributions to matrix A from the conservation of mass equations
    CM_eq_num = range(21,31)
    Q_left_i = [13, 14, 15, 20, 21, 22, 16, 17, 18, 19] # all +1 coefficients
    Q_right_i = [15, 17, 19, 22, 24, 26, 27, 28, 29, 30] # all -1 coefficients
    Q_mid_i = [14, 16, 18, 21, 23, 25, 23, 24, 25, 26] # first 6 elements -1 coefficients, last four +1
    for i in range(0,len(CM_eq_num)):
        M_Q[CM_eq_num[i], Q_left_i[i]] = 1
        M_Q[CM_eq_num[i], Q_right_i[i]] = -1
        if CM_eq_num[i]<=26:
            M_Q[CM_eq_num[i], Q_mid_i[i]] = -1
        else:
            M_Q[CM_eq_num[i], Q_mid_i[i]] = 1

    #   Add pressure and flow rate matrix contributions
    A = M_P + M_Q
    return A

def get_reduced_system(R, P1, P2, rho, dh):

    g = 9.81  # m/s^2
    rho1 = rho[0,0]  # kg/m^3
    rho2 = rho[0,1]  # kg/m^3
    rho3 = rho[0,2]  # kg/m^3

    #   Change in heights
    dh_in = dh[0,0]  # m
    dh_out = dh[0,1]  # m

    #   Create A matrix and b vector
    A = np.zeros([7,7])
    b = np.zeros([7,1])

    #   Define coefficients
    A[0,0] = 1
    A[1,1] = 1
    A[2,3] = 1
    A[3,4] = 1
    A[3,5] = rho[0,1]/rho[0,0]
    A[3,6] = -rho[0,2]/rho[0,0]
    A[4,0] = 1
    A[4,2] = -1
    A[4,4] = -R[0,0]
    A[5,1] = 1
    A[5,2] = -1
    A[5,5] = -R[0,1]
    A[6,2] = 1
    A[6,3] = -1
    A[6,6] = -R[0,2]

    #   Define b vector values
    b[0,0] = P1
    b[1,0] = P2
    b[4,0] = -rho1*g*dh_in
    b[5,0] = -rho2*g*dh_in
    b[6,0] = -rho3*g*dh_out

    return A, b

def get_reduced_output(A, b):

    #   Solve system
    x = np.linalg.solve(A, b)
    print(x)
    print('Q_out = ', (6*10**7)*x[6,0]/4)

def get_full_b_vector(P1, P2):
    """Define the b vector representing the RHS of the linear system."""

    b = np.zeros([31, 1])
    b[0, 0] = P1
    b[1, 0] = P2
    return b

def get_full_output(A, b):
    """Solve the linear system Ax=b and separate the output array x into pressure and flow rate arrays."""

    #   Solve linear system
    x = np.linalg.solve(A, b)
    P = x[0:13,0]
    Q = x[13:,0]

    #   Convert units
    P = P/100       #   Pa -> mbar
    Q = (6*10**7)*Q  #   m^3/s -> mL/min

    return P, Q

def print_pressures(pressures, pressure_labels):
    pressures = np.asarray(pressures).flatten()

    if len(pressures) != 13:
        raise ValueError("The pressure array must contain exactly 13 values.")

    if len(pressure_labels) != 13:
        raise ValueError("You must provide exactly 13 pressure labels.")

    segments = [
        ("Sol. 1 Branch", pressures[0:4], pressure_labels[0:4]),
        ("Sol. 2 Branch", pressures[4:8], pressure_labels[4:8]),
        ("Output", pressures[8:13], pressure_labels[8:13]),
    ]

    label_width = 18
    value_width = 10
    column_width = label_width + value_width + 1
    gap = "    "
    number_of_columns = len(segments)

    total_width = (
            number_of_columns * column_width
            + (number_of_columns - 1) * len(gap)
    )

    section_break = "-" * total_width

    print(section_break)
    # Segment headings
    print(gap.join(
        name.center(column_width)
        for name, _, _ in segments
    ))

    # Separator lines
    print(gap.join("-" * column_width for _ in segments))

    # Column headings
    heading = f"{'Location':<{label_width}} {'Pressure':>{value_width}}"
    print(gap.join(heading for _ in segments))

    # Unit labels
    units = f"{'':<{label_width}} {'(mbar)':>{value_width}}"
    print(gap.join(units for _ in segments))

    # Separator lines
    print(gap.join("-" * column_width for _ in segments))

    maximum_rows = max(len(values) for _, values, _ in segments)

    for row in range(maximum_rows):
        formatted_columns = []

        for _, values, labels in segments:
            if row < len(values):
                entry = (
                    f"{labels[row]:<{label_width}} "
                    f"{values[row]:>{value_width}.1f}"
                )
            else:
                entry = " " * column_width

            formatted_columns.append(entry)

        print(gap.join(formatted_columns))

    print(section_break)

def print_flow_rates(flow_rates, flow_labels):

    total_flow_rate = sum(flow_rates[14:])
    flow_unit = "mL/min"

    flow_rates = np.asarray(flow_rates).flatten()

    if len(flow_rates) != 18:
        raise ValueError(
            "The flow-rate array must contain exactly 18 values."
        )

    if len(flow_labels) != 18:
        raise ValueError(
            "You must provide exactly 18 flow-rate labels."
        )

    segments = [
        (
            "Sol. 1 Branch",
            flow_rates[0:7],
            flow_labels[0:7]
        ),
        (
            "Sol. 2 Branch",
            flow_rates[7:14],
            flow_labels[7:14]
        ),
        (
            "Output",
            flow_rates[14:18],
            flow_labels[14:18]
        ),
    ]

    label_width = 18
    value_width = 10
    column_width = label_width + value_width + 1
    gap = "    "

    number_of_columns = len(segments)

    total_width = (
        number_of_columns * column_width
        + (number_of_columns - 1) * len(gap)
    )

    section_break = "-" * total_width

    # Top section break
    print(section_break)

    # Segment headings
    print(gap.join(
        name.center(column_width)
        for name, _, _ in segments
    ))

    # Separator lines
    print(gap.join("-" * column_width for _ in segments))

    # Column headings
    heading = (
        f"{'Resistor':<{label_width}} "
        f"{'Flow Rate':>{value_width}}"
    )
    print(gap.join(heading for _ in segments))

    # Units
    units = (
        f"{'':<{label_width}} "
        f"{f'({flow_unit})':>{value_width}}"
    )
    print(gap.join(units for _ in segments))

    # Separator lines
    print(gap.join("-" * column_width for _ in segments))

    # Values
    maximum_rows = max(
        len(values) for _, values, _ in segments
    )

    for row in range(maximum_rows):
        formatted_columns = []

        for _, values, labels in segments:
            if row < len(values):
                entry = (
                    f"{labels[row]:<{label_width}} "
                    f"{values[row]:>{value_width}.3f}"
                )
            else:
                entry = " " * column_width

            formatted_columns.append(entry)

        print(gap.join(formatted_columns))

    # Bottom section break
    print(section_break)

    # Total flow rate
    print(
        f"{'Total Flow Rate:':<{label_width}} "
        f"{total_flow_rate:.3f} {flow_unit}"
    )

def print_full_result(p_array, q_array):
    """Print the output of the circuit analysis."""

    #   Create data labels
    p_labels = ['Sol. 1 Inlet', 'Junct. 1', 'Junct. 2', 'Junct. 3',
                'Sol. 2 Inlet', 'Junct. 1', 'Junct. 2', 'Junct. 3',
                'Junct. 1', 'Junct. 2', 'Junct. 3', 'Junct. 4', 'Outlet']
    q_labels = ['1', '2', '3', '4', '5', '6', '7', '1', '2', '3', '4', '5', '6', '7','1', '2', '3', '4']

    print('='*95)
    print('FULL RESULTS')
    print_pressures(p_array, p_labels)
    print_flow_rates(q_array, q_labels)
    print('\n')
    print('=' * 95)

    #   Estimate change in radius due to pressure
    E = 4 * (10 ** 6)
    nu = 0.44
    r_i = (794 / 2) * 10 ** -6
    r_o = (2.38125 / 2) * 10 ** -3
    dr = lambda p: ((1 - nu) * ((p * r_i ** 2) / (r_o ** 2 - r_i ** 2)) * r_i + ((1 + nu) / r_i) * (r_i ** 2) * (
                r_o ** 2) * p / (r_o ** 2 - r_i ** 2)) / E

    dr_est = dr(100 * 100)
    print(dr_est * 10 ** 6)

    # A = get_full_matrix(R)
    # b = get_full_b_vector(P01, P02)
    # P, Q = get_full_output(A,b)
    # print_full_result(P, Q)
    # get_reduced_system(R, P01, P02)

    D = 0.794 * 10 ** -3  # m
    v = (Qout / 4) / (np.pi * (D / 2) ** 2)  # m/s
    Re = rho[0, 2] * v * D / 0.001002

    print(Re)