"""Plot per-outlet HTR flow rate over one or two inlet-pressure ranges.

Pressure values supplied to this script are in mbar.  They are converted to
Pa before being passed to HTR_module_circuit_analysis.main().  Flow rates
returned by that function are in mL/min.
"""

from __future__ import annotations

from collections.abc import Sequence
import tkinter as tk
from tkinter import messagebox, ttk

import matplotlib.pyplot as plt
import numpy as np

from HTR_module_circuit_analysis import main as calculate_htr_flow


MBAR_TO_PA = 100.0
INVALID_FLOW_RATE = -1.0


def _as_pressure_array(pressure_range_mbar: Sequence[float]) -> np.ndarray:
    """Return a validated, one-dimensional pressure array."""
    pressures = np.asarray(pressure_range_mbar, dtype=float)
    if pressures.ndim != 1 or pressures.size == 0:
        raise ValueError("A pressure range must be a non-empty 1-D sequence.")
    if not np.all(np.isfinite(pressures)):
        raise ValueError("Pressure ranges may only contain finite values.")
    return pressures


def _per_outlet_flow_rate(p1_mbar: float, p2_mbar: float) -> float:
    """Calculate the flow through one of the four outlets in mL/min."""
    q = calculate_htr_flow(p1_mbar * MBAR_TO_PA, p2_mbar * MBAR_TO_PA)
    flow_rates = np.asarray(q, dtype=float).ravel()
    if flow_rates.size < 3:
        raise ValueError("Q must contain at least three flow-rate values.")
    if abs(flow_rates[0]) + abs(flow_rates[1]) > flow_rates[-1]:
        return INVALID_FLOW_RATE
    return float(flow_rates[-1]) / 4.0


def calculate_case_1(pressure_range_mbar: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Calculate flow when both inlet pressures are equal."""
    pressures = _as_pressure_array(pressure_range_mbar)
    flow_rates = np.array(
        [_per_outlet_flow_rate(pressure, pressure) for pressure in pressures]
    )
    return pressures, flow_rates


def calculate_case_2(
    pressure_1_range_mbar: Sequence[float],
    pressure_2_range_mbar: Sequence[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate flow for every combination of the two inlet-pressure ranges."""
    pressure_1 = _as_pressure_array(pressure_1_range_mbar)
    pressure_2 = _as_pressure_array(pressure_2_range_mbar)

    flow_rates = np.empty((pressure_2.size, pressure_1.size), dtype=float)
    for row, p2 in enumerate(pressure_2):
        for column, p1 in enumerate(pressure_1):
            flow_rates[row, column] = _per_outlet_flow_rate(p1, p2)
    return pressure_1, pressure_2, flow_rates


def plot_case_1(pressure_range_mbar: Sequence[float]) -> tuple[plt.Figure, plt.Axes]:
    """Plot per-outlet flow rate for equal inlet pressures."""
    pressures, flow_rates = calculate_case_1(pressure_range_mbar)
    figure, axes = plt.subplots()
    axes.plot(pressures, flow_rates, marker="o")
    axes.set(
        xlabel="Both inlet pressures (mbar)",
        ylabel="Flow rate per outlet (mL/min)",
        title="Case 1: Equal inlet pressures",
    )
    axes.grid(True, alpha=0.3)
    figure.tight_layout()
    return figure, axes


def plot_case_2(
    pressure_1_range_mbar: Sequence[float],
    pressure_2_range_mbar: Sequence[float],
) -> tuple[plt.Figure, plt.Axes]:
    """Plot a filled contour of all case-2 pressure combinations."""
    pressure_1, pressure_2, flow_rates = calculate_case_2(
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
        title="Case 2: Flow rate for all inlet-pressure combinations",
    )
    figure.tight_layout()
    return figure, axes


def _pressure_range(start: float, stop: float, points: int) -> np.ndarray:
    if points < 2:
        raise ValueError("The number of points must be at least 2.")
    return np.linspace(start, stop, points)


def prompt_for_inputs() -> tuple[int, np.ndarray, np.ndarray | None] | None:
    """Show an input form and return the selected case and pressure ranges."""
    root = tk.Tk()
    root.title("HTR Flow Rate Plot")
    root.resizable(False, False)

    case = tk.IntVar(value=1)
    p1_start = tk.StringVar(value="100")
    p1_stop = tk.StringVar(value="500")
    p2_start = tk.StringVar(value="100")
    p2_stop = tk.StringVar(value="500")
    points = tk.StringVar(value="25")
    result: dict[str, object] = {}

    frame = ttk.Frame(root, padding=16)
    frame.grid(row=0, column=0, sticky="nsew")
    ttk.Label(frame, text="Select a calculation case:").grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 4)
    )
    ttk.Radiobutton(
        frame, text="Case 1: equal inlet pressures", variable=case, value=1
    ).grid(row=1, column=0, columnspan=3, sticky="w")
    ttk.Radiobutton(
        frame, text="Case 2: different inlet pressures", variable=case, value=2
    ).grid(row=2, column=0, columnspan=3, sticky="w", pady=(0, 12))

    ttk.Label(frame, text="Pressure range").grid(row=3, column=0, sticky="w")
    ttk.Label(frame, text="Start (mbar)").grid(row=3, column=1, padx=6)
    ttk.Label(frame, text="Stop (mbar)").grid(row=3, column=2, padx=6)

    ttk.Label(frame, text="Inlet pressure 1").grid(row=4, column=0, sticky="w")
    ttk.Entry(frame, textvariable=p1_start, width=14).grid(row=4, column=1, padx=6, pady=3)
    ttk.Entry(frame, textvariable=p1_stop, width=14).grid(row=4, column=2, padx=6, pady=3)

    ttk.Label(frame, text="Inlet pressure 2").grid(row=5, column=0, sticky="w")
    p2_start_entry = ttk.Entry(frame, textvariable=p2_start, width=14)
    p2_start_entry.grid(row=5, column=1, padx=6, pady=3)
    p2_stop_entry = ttk.Entry(frame, textvariable=p2_stop, width=14)
    p2_stop_entry.grid(row=5, column=2, padx=6, pady=3)

    ttk.Label(frame, text="Points per pressure range").grid(
        row=6, column=0, sticky="w", pady=(8, 0)
    )
    ttk.Entry(frame, textvariable=points, width=14).grid(
        row=6, column=1, padx=6, pady=(8, 0)
    )

    def update_case_fields(*_args: object) -> None:
        state = "normal" if case.get() == 2 else "disabled"
        p2_start_entry.configure(state=state)
        p2_stop_entry.configure(state=state)

    def submit() -> None:
        try:
            point_count = int(points.get())
            pressure_1 = _pressure_range(
                float(p1_start.get()), float(p1_stop.get()), point_count
            )
            pressure_2 = None
            if case.get() == 2:
                pressure_2 = _pressure_range(
                    float(p2_start.get()), float(p2_stop.get()), point_count
                )
        except ValueError as error:
            messagebox.showerror(
                "Invalid input",
                f"Enter numeric pressure limits and a whole-number point count.\n\n{error}",
                parent=root,
            )
            return

        result["values"] = (case.get(), pressure_1, pressure_2)
        root.destroy()

    def cancel() -> None:
        root.destroy()

    button_frame = ttk.Frame(frame)
    button_frame.grid(row=7, column=0, columnspan=3, sticky="e", pady=(16, 0))
    ttk.Button(button_frame, text="Cancel", command=cancel).grid(row=0, column=0, padx=4)
    ttk.Button(button_frame, text="Create Plot", command=submit).grid(row=0, column=1)

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
    return result.get("values")  # type: ignore[return-value]


def main() -> None:
    inputs = prompt_for_inputs()
    if inputs is None:
        return

    selected_case, pressure_1, pressure_2 = inputs
    if selected_case == 1:
        figure, _ = plot_case_1(pressure_1)
    else:
        if pressure_2 is None:
            raise RuntimeError("Case 2 requires a second pressure range.")
        figure, _ = plot_case_2(pressure_1, pressure_2)

    plt.show()


if __name__ == "__main__":
    main()
