"""Fit tab-separated fluid property tables with polynomials in temperature.

Run this script to select a file, property, and polynomial degree(s) in a GUI.
Case 1 overlays any comma-separated list of degrees on the measured data.
Case 2 plots one fit and prints '"property name": "equation"' to the console and GUI.
Both cases display their plots inside the GUI window.
T and the fitted property retain the units in the input file's headers.
Requires NumPy and Matplotlib; Tkinter is included with standard Windows Python.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tkinter as tk
from tkinter import messagebox, ttk
import warnings

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from numpy.polynomial import Polynomial


DATA_DIRECTORY = Path(__file__).resolve().parent / "Fluid Property Data"


def read_property_table(path: Path) -> tuple[str, dict[str, np.ndarray]]:
    """Read numeric columns, preserving headers and rejecting malformed rows."""
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        rows = [row for row in reader if any(cell.strip() for cell in row)]
    if len(rows) < 2:
        raise ValueError("The file must contain a header and at least one data row.")
    headers = [cell.strip() for cell in rows[0]]
    if any(not name for name in headers) or len(set(headers)) != len(headers):
        raise ValueError("Column headers must be nonempty and unique.")
    if any(len(row) != len(headers) for row in rows[1:]):
        raise ValueError("Each data row must have the same number of tabs as the header.")
    columns = {}
    for index, name in enumerate(headers):
        try:
            values = np.array([float(row[index]) for row in rows[1:]])
        except ValueError:
            continue  # Text columns, such as Phase, are not fit candidates.
        if np.all(np.isfinite(values)):
            columns[name] = values
    temperatures = [name for name in columns if name.lower().startswith("temperature")]
    if len(temperatures) != 1:
        raise ValueError("Expected one finite numeric column named Temperature (with units).")
    if len(columns) < 2:
        raise ValueError("No finite numeric property columns were found.")
    return temperatures[0], columns


def parse_orders(text: str, compare: bool, max_order: int) -> list[int]:
    """Validate polynomial degrees (degree zero is a constant fit)."""
    try:
        orders = [int(part.strip()) for part in text.split(",")]
    except ValueError as error:
        raise ValueError("Enter whole-number orders, e.g. 1, 2, 3 or a single 3.") from error
    if not compare and len(orders) != 1:
        raise ValueError("Case 2 requires exactly one order.")
    if any(order < 0 or order > max_order for order in orders):
        raise ValueError(f"Orders must be between 0 and {max_order} for these temperatures.")
    return list(dict.fromkeys(orders))


def fit_polynomial(temperature: np.ndarray, values: np.ndarray, order: int) -> Polynomial:
    """Scale temperature during least squares to improve numerical conditioning."""
    if order < 0 or order >= np.unique(temperature).size:
        raise ValueError("The order must be less than the number of unique temperatures.")
    if order == 0:
        return Polynomial([float(np.mean(values))])
    with warnings.catch_warnings():
        warnings.simplefilter("error", np.exceptions.RankWarning)
        try:
            model = Polynomial.fit(temperature, values, order)
        except np.exceptions.RankWarning as error:
            raise ValueError("This fit is poorly conditioned. Try a lower order.") from error
    if not np.all(np.isfinite(model.coef)):
        raise ValueError("The fit produced nonfinite coefficients. Try a lower order.")
    return model


def regression_equation(model: Polynomial) -> str:
    """Return a Python-compatible polynomial in the original temperature units."""
    coefficients = model.convert().coef
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("Cannot express this fit in powers of T. Try a lower order.")
    terms = []
    for power in range(len(coefficients) - 1, -1, -1):
        coefficient = coefficients[power]
        if coefficient == 0:
            continue
        term = f"{abs(coefficient):.16g}"
        if power == 1:
            term += "*T"
        elif power > 1:
            term += f"*T**{power}"
        sign = "-" if coefficient < 0 else "+"
        terms.append((sign, term))
    if not terms:
        return "0"
    sign, term = terms[0]
    return ("-" if sign == "-" else "") + term + "".join(
        f" {sign} {term}" for sign, term in terms[1:]
    )


def plot_comparison(temperature, values, models, temperature_name, property_name, filename,
                    figure=None):
    """Overlay each requested fit and label it with its training-data RMSE."""
    if figure is None:
        figure = Figure(figsize=(9, 5))
    figure.clear()
    axis = figure.add_subplot(111)
    axis.scatter(temperature, values, color="black", s=18, label="Source data", zorder=3)
    grid = np.linspace(temperature.min(), temperature.max(), 500)
    for order, model in models.items():
        rmse = np.sqrt(np.mean((model(temperature) - values) ** 2))
        axis.plot(grid, model(grid), label=f"Order {order}; RMSE = {rmse:.4g}")
    axis.set(xlabel=temperature_name, ylabel=property_name,
             title=f"{filename}\n{property_name} vs. temperature")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.tight_layout()
    return figure


def main() -> None:
    root = tk.Tk()
    root.title("Temperature-dependent property fitting")
    frame = ttk.Frame(root, padding=16)
    frame.grid(sticky="nsew")
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    frame.columnconfigure(1, weight=1)
    frame.rowconfigure(9, weight=1)

    files = sorted(DATA_DIRECTORY.glob("*.txt"), key=lambda path: path.name.lower())
    if not files:
        messagebox.showerror("No data files", f"No .txt files found in:\n{DATA_DIRECTORY}", parent=root)
        root.destroy()
        return

    file_name = tk.StringVar(value=files[0].name)
    property_name = tk.StringVar()
    mode = tk.IntVar(value=1)
    orders = tk.StringVar(value="1, 2, 3")
    info = tk.StringVar()
    table = {}
    temperature_name = ""

    ttk.Label(frame, text="Data file:").grid(row=0, column=0, sticky="w")
    file_box = ttk.Combobox(frame, textvariable=file_name, state="readonly",
                            values=[path.name for path in files], width=50)
    file_box.grid(row=0, column=1, sticky="ew", pady=4)
    ttk.Label(frame, text="Property:").grid(row=1, column=0, sticky="w")
    property_box = ttk.Combobox(frame, textvariable=property_name, state="readonly")
    property_box.grid(row=1, column=1, sticky="ew", pady=4)
    ttk.Label(frame, textvariable=info, wraplength=620).grid(
        row=2, column=0, columnspan=2, sticky="w", pady=8)

    def update_mode():
        orders.set("1, 2, 3" if mode.get() == 1 else "3")

    ttk.Radiobutton(frame, text="Case 1: Compare multiple orders on a plot",
                    variable=mode, value=1, command=update_mode).grid(
        row=3, column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(frame, text="Case 2: Plot and print the equation for one order",
                    variable=mode, value=2, command=update_mode).grid(
        row=4, column=0, columnspan=2, sticky="w")
    ttk.Label(frame, text="Polynomial order(s):").grid(row=5, column=0, sticky="w")
    ttk.Entry(frame, textvariable=orders).grid(row=5, column=1, sticky="ew", pady=8)
    ttk.Label(frame, text="Case 1: any comma-separated list (e.g. 1, 2, 3, 5). Case 2: one integer.").grid(
        row=6, column=0, columnspan=2, sticky="w")
    output = tk.Text(frame, width=85, height=4, wrap="word")
    output.grid(row=8, column=0, columnspan=2, sticky="ew", pady=10)

    plot_frame = ttk.Frame(frame)
    plot_frame.grid(row=9, column=0, columnspan=2, sticky="nsew")
    plot_frame.columnconfigure(0, weight=1)
    plot_frame.rowconfigure(0, weight=1)
    figure = Figure(figsize=(9, 4))
    canvas = FigureCanvasTkAgg(figure, master=plot_frame)
    canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
    toolbar = NavigationToolbar2Tk(canvas, plot_frame, pack_toolbar=False)
    toolbar.grid(row=1, column=0, sticky="ew")

    def load_file(_event=None):
        nonlocal table, temperature_name
        table = {}
        property_box["values"] = ()
        property_name.set("")
        output.delete("1.0", tk.END)
        figure.clear()
        toolbar.update()
        canvas.draw()
        try:
            temperature_name, table = read_property_table(DATA_DIRECTORY / file_name.get())
        except (OSError, ValueError) as error:
            info.set("Unable to load this file.")
            messagebox.showerror("Invalid file", str(error), parent=root)
            return
        properties = [name for name in table if name != temperature_name]
        property_box["values"] = properties
        property_name.set(properties[0])
        t = table[temperature_name]
        info.set(f"T uses {temperature_name}. Data range: {t.min():g} to {t.max():g}; "
                 f"{t.size} rows. Maximum order: {np.unique(t).size - 1}.")

    def run_fit():
        try:
            if not table or property_name.get() not in table:
                raise ValueError("Select a valid file and numeric property first.")
            t = table[temperature_name]
            values = table[property_name.get()]
            degrees = parse_orders(orders.get(), mode.get() == 1, np.unique(t).size - 1)
            models = {degree: fit_polynomial(t, values, degree) for degree in degrees}
            output.delete("1.0", tk.END)
            if mode.get() == 1:
                output.insert(tk.END, "Compare the curves, then select Case 2 and enter your chosen order.")
            else:
                equation = regression_equation(models[degrees[0]])
                line = f"{json.dumps(property_name.get())}: {json.dumps(equation)}"
                print(line, flush=True)
                output.insert(tk.END, line + f"\n\nT: {temperature_name}; fitted range: {t.min():g} to {t.max():g}.")
            plot_comparison(t, values, models, temperature_name, property_name.get(),
                            file_name.get(), figure=figure)
            toolbar.update()
            canvas.draw()
        except (ValueError, np.linalg.LinAlgError) as error:
            messagebox.showerror("Cannot fit property", str(error), parent=root)

    ttk.Button(frame, text="Fit property", command=run_fit).grid(row=7, column=1, sticky="e", pady=8)
    file_box.bind("<<ComboboxSelected>>", load_file)

    def close():
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", close)
    load_file()
    root.mainloop()


if __name__ == "__main__":
    main()
