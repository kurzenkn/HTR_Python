"""Calculate and display HTR results from user-supplied inlet pressures.

This standalone script combines the circuit model and pressure-range plotting
interface. Workbook parameters are loaded once per analysis, before any loop
over pressure values.
"""

from __future__ import annotations

from collections.abc import Sequence
import tkinter as tk
from tkinter import messagebox, ttk

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MBAR_TO_PA = 100.0
INVALID_FLOW_RATE = -1.0
INVALID_FLOW_RATIO = -1.0


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


def _as_pressure_array(pressure_range_mbar: Sequence[float]) -> np.ndarray:
    """Return a validated, one-dimensional pressure array."""
    pressures = np.asarray(pressure_range_mbar, dtype=float)
    if pressures.ndim != 1 or pressures.size == 0:
        raise ValueError("A pressure range must be a non-empty 1-D sequence.")
    if not np.all(np.isfinite(pressures)):
        raise ValueError("Pressure ranges may only contain finite values.")
    return pressures


def _as_pressure_value(pressure_mbar: float) -> float:
    """Return a validated scalar pressure value."""
    pressure = float(pressure_mbar)
    if not np.isfinite(pressure):
        raise ValueError("Pressure values must be finite.")
    return pressure


def _load_circuit_parameters() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read workbook values once for a complete analysis."""
    resistances = calc_resistors()
    reduced_resistances = get_reduced_resistances(resistances)
    densities = get_densities()
    heights = get_heights()
    return reduced_resistances, densities, heights


def calculate_single_case(
    p1_mbar: float,
    p2_mbar: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate pressure and flow results for one inlet-pressure pair."""
    pressure_1 = _as_pressure_value(p1_mbar)
    pressure_2 = _as_pressure_value(p2_mbar)
    reduced_resistances, densities, heights = _load_circuit_parameters()
    flow_rates, pressures = calc_reduced_result(
        reduced_resistances,
        densities,
        pressure_1 * MBAR_TO_PA,
        pressure_2 * MBAR_TO_PA,
        heights,
    )
    return flow_rates, pressures, reduced_resistances


def _per_outlet_flow_rate(
    p1_mbar: float,
    p2_mbar: float,
    circuit_parameters: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> float:
    """Calculate the flow through one of the four outlets in mL/min."""
    reduced_resistances, densities, heights = circuit_parameters
    q, _ = calc_reduced_result(
        reduced_resistances,
        densities,
        p1_mbar * MBAR_TO_PA,
        p2_mbar * MBAR_TO_PA,
        heights,
    )
    flow_rates = np.asarray(q, dtype=float).ravel()
    if flow_rates.size < 3:
        raise ValueError("Q must contain at least three flow-rate values.")
    if abs(flow_rates[0]) + abs(flow_rates[1]) > flow_rates[-1]:
        return INVALID_FLOW_RATE
    return float(flow_rates[-1]) / 4.0


def calculate_case_2(pressure_range_mbar: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Calculate flow when both inlet pressures follow the same range."""
    pressures = _as_pressure_array(pressure_range_mbar)
    circuit_parameters = _load_circuit_parameters()
    flow_rates = np.array(
        [
            _per_outlet_flow_rate(pressure, pressure, circuit_parameters)
            for pressure in pressures
        ]
    )
    return pressures, flow_rates


def calculate_case_3(
    pressure_1_range_mbar: Sequence[float],
    pressure_2_range_mbar: Sequence[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return both pressure axes, per-outlet flows, and Q2/Q1 ratios."""
    pressure_1 = _as_pressure_array(pressure_1_range_mbar)
    pressure_2 = _as_pressure_array(pressure_2_range_mbar)
    reduced_resistances, densities, heights = _load_circuit_parameters()

    flow_rates = np.empty((pressure_2.size, pressure_1.size), dtype=float)
    flow_ratios = np.full_like(flow_rates, np.nan)
    for row, p2 in enumerate(pressure_2):
        for column, p1 in enumerate(pressure_1):
            q, _ = calc_reduced_result(
                reduced_resistances, densities,
                p1 * MBAR_TO_PA, p2 * MBAR_TO_PA, heights
            )
            q1, q2, qout = np.asarray(q, dtype=float).ravel()
            flow_rates[row, column] = (
                INVALID_FLOW_RATE if abs(q1) + abs(q2) > qout else qout / 4.0
            )
            # Replace negative ratios with the dummy value; zero Q1 is undefined.
            if q1 != 0:
                with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                    ratio = q2 / q1
                flow_ratios[row, column] = INVALID_FLOW_RATIO if ratio < 0 else ratio
    return pressure_1, pressure_2, flow_rates, flow_ratios


def plot_case_2(pressure_range_mbar: Sequence[float]) -> tuple[plt.Figure, plt.Axes]:
    """Plot per-outlet flow rate for equal inlet pressures."""
    pressures, flow_rates = calculate_case_2(pressure_range_mbar)
    figure, axes = plt.subplots()
    axes.plot(pressures, flow_rates, marker="o")
    axes.set(
        xlabel="Both inlet pressures (mbar)",
        ylabel="Flow rate per outlet (mL/min)",
        title="Case 2: Equal inlet pressures",
    )
    axes.grid(True, alpha=0.3)
    figure.tight_layout()
    return figure, axes


def plot_case_3(
    pressure_1_range_mbar: Sequence[float],
    pressure_2_range_mbar: Sequence[float],
) -> tuple[plt.Figure, plt.Axes]:
    """Show flow and Q2/Q1 contours; return the original flow figure and axes."""
    pressure_1, pressure_2, flow_rates, flow_ratios = calculate_case_3(
        pressure_1_range_mbar, pressure_2_range_mbar
    )
    pressure_1_grid, pressure_2_grid = np.meshgrid(pressure_1, pressure_2)

    figure, axes = plt.subplots(figsize=(8, 6))
    contour = axes.contourf(
        pressure_1_grid, pressure_2_grid, flow_rates, levels=20, cmap="viridis"
    )
    colorbar = figure.colorbar(contour, ax=axes)
    colorbar.set_label("Flow rate per outlet (mL/min)")
    axes.set(
        xlabel="Inlet pressure 1 (mbar)",
        ylabel="Inlet pressure 2 (mbar)",
        title="Case 3: Flow rate for all inlet-pressure combinations",
    )
    figure.tight_layout()

    ratio_figure, ratio_axes = plt.subplots(figsize=(8, 6))
    masked_ratios = np.ma.masked_invalid(flow_ratios)
    if masked_ratios.count():
        ratio_contour = ratio_axes.contourf(
            pressure_1_grid, pressure_2_grid, masked_ratios,
            levels=20, cmap="viridis",
        )
        ratio_colorbar = ratio_figure.colorbar(ratio_contour, ax=ratio_axes)
        ratio_colorbar.set_label("Flow rate ratio Q2 / Q1")
    else:
        ratio_axes.text(
            0.5, 0.5, "Q2 / Q1 is undefined for all selected pressures.",
            ha="center", va="center", transform=ratio_axes.transAxes,
        )
    ratio_axes.set(
        xlabel="Inlet pressure 1 (mbar)",
        ylabel="Inlet pressure 2 (mbar)",
        title="Case 3: Flow rate ratio Q2 / Q1",
    )
    ratio_figure.tight_layout()
    return figure, axes


def _pressure_range(start: float, stop: float, points: int) -> np.ndarray:
    if points < 2:
        raise ValueError("The number of points must be at least 2.")
    if not np.isfinite(start) or not np.isfinite(stop):
        raise ValueError("Pressure limits must be finite.")
    return np.linspace(start, stop, points)


InputValue = float | np.ndarray
InputSelection = tuple[int, InputValue, InputValue | None]


def prompt_for_inputs() -> InputSelection | None:
    """Show an input form and return the selected case and pressure inputs."""
    root = tk.Tk()
    root.title("HTR Analysis")
    root.resizable(False, False)

    case = tk.IntVar(value=1)
    p1_start = tk.StringVar(value="10")
    p1_stop = tk.StringVar(value="50")
    p2_start = tk.StringVar(value="10")
    p2_stop = tk.StringVar(value="50")
    points = tk.StringVar(value="25")
    result: dict[str, InputSelection] = {}

    frame = ttk.Frame(root, padding=16)
    frame.grid(row=0, column=0, sticky="nsew")
    ttk.Label(frame, text="Select an analysis case:").grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 4)
    )
    ttk.Radiobutton(
        frame,
        text="Case 1: one inlet-pressure pair (print results)",
        variable=case,
        value=1,
    ).grid(row=1, column=0, columnspan=3, sticky="w")
    ttk.Radiobutton(
        frame,
        text="Case 2: equal inlet-pressure range (plot)",
        variable=case,
        value=2,
    ).grid(row=2, column=0, columnspan=3, sticky="w")
    ttk.Radiobutton(
        frame,
        text="Case 3: different inlet-pressure ranges (plot)",
        variable=case,
        value=3,
    ).grid(row=3, column=0, columnspan=3, sticky="w", pady=(0, 12))

    ttk.Label(frame, text="Pressure input").grid(row=4, column=0, sticky="w")
    ttk.Label(frame, text="Value / start (mbar)").grid(row=4, column=1, padx=6)
    ttk.Label(frame, text="Stop (mbar)").grid(row=4, column=2, padx=6)

    ttk.Label(frame, text="Inlet pressure 1").grid(row=5, column=0, sticky="w")
    p1_start_entry = ttk.Entry(frame, textvariable=p1_start, width=18)
    p1_start_entry.grid(row=5, column=1, padx=6, pady=3)
    p1_stop_entry = ttk.Entry(frame, textvariable=p1_stop, width=14)
    p1_stop_entry.grid(row=5, column=2, padx=6, pady=3)

    ttk.Label(frame, text="Inlet pressure 2").grid(row=6, column=0, sticky="w")
    p2_start_entry = ttk.Entry(frame, textvariable=p2_start, width=18)
    p2_start_entry.grid(row=6, column=1, padx=6, pady=3)
    p2_stop_entry = ttk.Entry(frame, textvariable=p2_stop, width=14)
    p2_stop_entry.grid(row=6, column=2, padx=6, pady=3)

    points_label = ttk.Label(frame, text="Points per pressure range")
    points_label.grid(row=7, column=0, sticky="w", pady=(8, 0))
    points_entry = ttk.Entry(frame, textvariable=points, width=18)
    points_entry.grid(row=7, column=1, padx=6, pady=(8, 0))

    def update_case_fields(*_args: object) -> None:
        selected_case = case.get()
        p1_stop_entry.configure(state="normal" if selected_case != 1 else "disabled")
        p2_start_entry.configure(state="normal" if selected_case in (1, 3) else "disabled")
        p2_stop_entry.configure(state="normal" if selected_case == 3 else "disabled")
        range_state = "normal" if selected_case in (2, 3) else "disabled"
        points_entry.configure(state=range_state)
        points_label.configure(state=range_state)

    def submit() -> None:
        try:
            selected_case = case.get()
            if selected_case == 1:
                pressure_1: InputValue = _as_pressure_value(float(p1_start.get()))
                pressure_2: InputValue | None = _as_pressure_value(float(p2_start.get()))
            else:
                point_count = int(points.get())
                pressure_1 = _pressure_range(
                    float(p1_start.get()), float(p1_stop.get()), point_count
                )
                pressure_2 = None
                if selected_case == 3:
                    pressure_2 = _pressure_range(
                        float(p2_start.get()), float(p2_stop.get()), point_count
                    )
        except ValueError as error:
            messagebox.showerror(
                "Invalid input",
                f"Enter valid numeric pressures and a whole-number point count "
                f"for range cases.\n\n{error}",
                parent=root,
            )
            return

        result["values"] = (selected_case, pressure_1, pressure_2)
        root.destroy()

    def cancel() -> None:
        root.destroy()

    button_frame = ttk.Frame(frame)
    button_frame.grid(row=8, column=0, columnspan=3, sticky="e", pady=(16, 0))
    ttk.Button(button_frame, text="Cancel", command=cancel).grid(row=0, column=0, padx=4)
    ttk.Button(button_frame, text="Run Analysis", command=submit).grid(row=0, column=1)

    case.trace_add("write", update_case_fields)
    update_case_fields()
    root.protocol("WM_DELETE_WINDOW", cancel)
    root.bind("<Return>", lambda _event: submit())
    root.bind("<Escape>", lambda _event: cancel())
    root.update_idletasks()
    root.geometry(
        f"+{(root.winfo_screenwidth() - root.winfo_reqwidth()) // 2}"
        f"+{(root.winfo_screenheight() - root.winfo_reqheight()) // 2}"
    )
    root.mainloop()
    return result.get("values")


def main() -> None:
    """Run the selected single-point or pressure-range analysis."""
    inputs = prompt_for_inputs()
    if inputs is None:
        return

    selected_case, pressure_1, pressure_2 = inputs
    if selected_case == 1:
        if not isinstance(pressure_1, float) or not isinstance(pressure_2, float):
            raise RuntimeError("Case 1 requires two individual pressure values.")
        flow_rates, pressures, resistances = calculate_single_case(
            pressure_1, pressure_2
        )
        print_reduced_results(pressures, flow_rates, resistances)
        return

    if not isinstance(pressure_1, np.ndarray):
        raise RuntimeError("Pressure-range cases require pressure arrays.")

    if selected_case == 2:
        plot_case_2(pressure_1)
    else:
        if not isinstance(pressure_2, np.ndarray):
            raise RuntimeError("Case 3 requires a second pressure range.")
        plot_case_3(pressure_1, pressure_2)

    plt.show()


if __name__ == "__main__":
    main()
