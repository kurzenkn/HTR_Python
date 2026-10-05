import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def calc_resistors():
    """Return array of resistor values calculated from input parameters from the HTR Microfluidics Excel file."""

    #   Get dynamic viscosity values
    eta_1 = pd.read_excel("HTR_parameters.xlsx",
                          "Resistor Parameters",
                          usecols='L', skiprows=2, nrows=1, header = None).iat[0,0]
    eta_2 = pd.read_excel("HTR_parameters.xlsx",
                          "Resistor Parameters",
                          usecols='L', skiprows=3, nrows=1, header = None).iat[0,0]
    eta_3 = pd.read_excel("HTR_parameters.xlsx",
                          "Resistor Parameters",
                          usecols='L', skiprows=4, nrows=1, header = None).iat[0,0]

    #   Get custom resitor values from Excel sheet
    Rc1 = pd.read_excel("HTR_parameters.xlsx",
                        "Resistor Parameters",
                        usecols='P', skiprows=2, nrows=1, header = None).iat[0,0]
    Rc2 = pd.read_excel("HTR_parameters.xlsx",
                        "Resistor Parameters",
                        usecols='P', skiprows=3, nrows=1, header=None).iat[0, 0]

    #   Convert units from GPa*s/m^3 to Pa*s/m^3
    Rc1 = (10**9)*Rc1
    Rc2 = (10**9)*Rc2

    #   Get tube length values
    L = pd.read_excel("HTR_parameters.xlsx","Resistor Parameters",
                      usecols='D', skiprows=2, nrows=18,header = None).to_numpy().T
    L = L/1000 # convert from mm to m

    #   Get tube radius values
    d = pd.read_excel("HTR_parameters.xlsx",
                        "Resistor Parameters",
                        usecols='E', skiprows=2, nrows=18, header=None).to_numpy().T
    r = d/(2*(10**6)) # convert diameter to radius and from microns to m

    #   Calculate tube cross-sectional area
    A = np.pi*r**2

    #   Calculate resistor values
    R = np.zeros([1, 18])
    R[0,0:7] = 8*eta_1*L[0,0:7]/(A[0,0:7]*r[0,0:7]**2)
    R[0,7:14] = 8*eta_2*L[0,7:14]/(A[0,7:14]*r[0,7:14]**2)
    R[0,14:] = 8*eta_3*L[0,14:]/(A[0,14:]*r[0,14:]**2)

    #   Add custom resistors in series with R11 and R12
    R[0,0] = R[0,0]+Rc1
    R[0,7] = R[0,7]+Rc2
    # print("R_h =", R)
    return R

def get_input_pressures():
    """Set the pressure values for the two reservoirs."""

    P1 = pd.read_excel("HTR_parameters.xlsx",
                          "Input",
                          usecols='C', skiprows=2, nrows=1, header=None).iat[0, 0]
    P2 = pd.read_excel("HTR_parameters.xlsx",
                          "Input",
                          usecols='C', skiprows=3, nrows=1, header=None).iat[0, 0]

    #   Convert input pressure from mbar to Pa
    P1 = 100*P1
    P2 = 100*P2

    return P1, P2

def get_reduced_resistances(R):
    """Define matrix A for the linear system of equations for the reduced model."""

    #   Calculate equivalent resistance values
    R_eq = lambda r1, r2, r3: r1 + (1/r2 + 1/r3)**-1

    R1_left = R_eq(R[0,1], R[0,3], R[0,4])
    R1_right = R_eq(R[0,2], R[0,5], R[0,6])
    R2_left = R_eq(R[0,8], R[0,10], R[0,11])
    R2_right = R_eq(R[0, 9], R[0, 12], R[0, 13])

    R1_eq = R_eq(R[0,0],R1_left, R1_right)
    R2_eq = R_eq(R[0,7], R2_left, R2_right)
    R_out = (1/R[0,14] + 1/R[0,15] + 1/R[0,16] + 1/R[0,17])**-1
    R_tot = R_eq(R_out, R1_eq, R2_eq)

    R_reduced = np.array([[R1_eq, R2_eq, R_out, R_tot]])

    return R_reduced

def get_densities():
    """Get density values for solution 1, solution 2, and combined solution from the Excel sheet."""

    rho = pd.read_excel("HTR_parameters.xlsx","Resistor Parameters",
                      usecols='L', skiprows=7, nrows=3,header = None).to_numpy().T

    return rho

def get_heights():
    dh = pd.read_excel("HTR_parameters.xlsx","Resistor Parameters",
                      usecols='P', skiprows=6, nrows=2,header = None).to_numpy().T
    dh = dh/1000    #   mm -> m
    return dh

def calc_reduced_result(R, rho, P1, P2, dh):
    """Calculate the flow rates for reduced system."""

    g = 9.81        #   m/s^2

    #   Define constants
    c1 = P1 + rho[0,0]*g*dh[0,0] +rho[0,2]*g*dh[0,1]
    a1 = R[0,0] + R[0,2]*rho[0,0]/rho[0,2]
    b1 = R[0,2]*rho[0,1]/rho[0,2]
    c2 = P2 + rho[0,1]*g*dh[0,0]+rho[0,2]*g*dh[0,1]
    a2 = R[0,2]*rho[0,0]/rho[0,2]
    b2 = R[0,1] + R[0,2]*rho[0,1]/rho[0,2]

    #   Calculate flow rates
    Q1 = (c1-c2*b1/b2)/(a1-a2*b1/b2)
    Q2 = (c2-a2*Q1)/b2
    Qout = (rho[0,0]*Q1 + rho[0,1]*Q2)/rho[0,2]

    Q = np.array([[Q1, Q2, Qout]])

    #   Calculate intermediate pressure
    p_inter = R[0,2]*Qout-rho[0,2]*g*dh[0,1]
    P = np.array([[P1, P2, p_inter, 0]])

    #   Convert units
    P = P/100  # Pa -> mbar
    Q = (6*10**7)*Q  # m^3/s -> mL/min

    return Q, P

def print_reduced_results(P,Q,R,):

    pressure_labels = ["Sol. 1 Inlet", "Sol. 2 Inlet", "T Junction", "Output"]
    flow_labels = ["Sol. 1 Branch", "Sol. 2 Branch", "Output", "Per Outlet"]
    pressure_unit = "mbar"
    flow_unit = "mL/min"

    R = R/(10**9)

    P = np.asarray(P).flatten()
    Q = np.asarray(Q).flatten()
    R = np.asarray(R).flatten()

    if len(P) != 4:
        raise ValueError("P must contain exactly 4 values.")

    if len(Q) != 3:
        raise ValueError("Q must contain exactly 3 values.")

    if len(R) != 4:
        raise ValueError("R must contain exactly 4 values.")

    if len(pressure_labels) != 4:
        raise ValueError("Four pressure labels are required.")



    flow_values = np.append(Q, Q[-1] / 4)

    resistor_labels = [
        "Sol. 1 Branch",
        "Sol. 2 Branch",
        "Output",
        "Total"
    ]

    segments = [
        (
            "Pressure Results",
            P,
            pressure_labels,
            "Location",
            "Pressure",
            pressure_unit
        ),
        (
            "Flow Rate Results",
            flow_values,
            flow_labels,
            "Location",
            "Flow Rate",
            flow_unit
        ),
        (
            "Resistance Results",
            R,
            resistor_labels,
            "Resistor",
            "Resistance",
            "GPa*s/m^3"
        ),
    ]

    label_width = 18
    value_width = 14
    column_width = label_width + value_width + 1
    gap = "    "

    total_width = (
        len(segments) * column_width
        + (len(segments) - 1) * len(gap)
    )

    section_break = "-" * total_width
    column_break = gap.join(
        "-" * column_width for _ in segments
    )

    # Top section break
    print(section_break)

    # Results titles
    print(gap.join(
        title.center(column_width)
        for title, _, _, _, _, _ in segments
    ))

    # Lines below results titles
    print(column_break)

    # Column headings
    headings = []

    for _, _, _, row_heading, value_heading, _ in segments:
        headings.append(
            f"{row_heading:<{label_width}} "
            f"{value_heading:>{value_width}}"
        )

    print(gap.join(headings))

    # Units
    units = []

    for _, _, _, _, _, unit in segments:
        units.append(
            f"{'':<{label_width}} "
            f"{f'({unit})':>{value_width}}"
        )

    print(gap.join(units))

    # Lines below column headings and units
    print(column_break)

    maximum_rows = max(
        len(values)
        for _, values, _, _, _, _ in segments
    )

    for row in range(maximum_rows):

        # Add a dashed line before the final row in every column
        if row == maximum_rows - 1:
            print(column_break)

        formatted_columns = []

        for _, values, labels, _, _, _ in segments:
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

def main(P_in1=None, P_in2=None):
    """Return calculated pressure and flow rate values."""
    R = calc_resistors()
    R_red = get_reduced_resistances(R)
    rho = get_densities()
    dh = get_heights()

    if P_in1 is None:
        P01, P02 = get_input_pressures()
        Q, P = calc_reduced_result(R_red, rho, P01, P02, dh)
        print_reduced_results(P, Q, R_red)
    elif P_in2 is None:
        P01 = P_in1
        P02 = P01
        Q, P = calc_reduced_result(R_red, rho, P01, P02, dh)
        return Q
    else:
        P01 = P_in1
        P02 = P_in2
        Q, P = calc_reduced_result(R_red, rho, P01, P02, dh)
        return Q


if __name__ == "__main__":
    main()