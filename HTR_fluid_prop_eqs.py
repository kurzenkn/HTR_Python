EQUATIONS = {
    "water_props" : {
        "name": "Water",
        "rho_func": lambda T: -8.659210819301207e-08*T**4 + 3.116127438234955e-05*T**3
                              - 0.006675870205665397*T**2 + 0.02446817448164573*T + 1000.149935458138,
        "eta_func": lambda T: 1.719022795589096e-11*T**4 - 5.401584091643156e-09*T**3
                             + 6.841381583626284e-07*T**2 - 4.554469055593262e-05*T + 0.001680424280620044,
        "inputs": {
            "T": {
                "min": 15,
                "max": 99,
                "units": "C",
            }
        },
        "rho_units": "kg/m^3",
        "eta_units": "Pa*s",
        "rho_eq_ref": "NIST Chemistry WebBook, fitted table values at 1 C increments.",
        "eta_eq_ref": "NIST Chemistry WebBook, fitted table values at 1 C increments.",
    }
}