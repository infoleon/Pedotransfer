# Pedotransfer Function

Pedotransfer function to predict field capacity effective saturation of an arable field from van Genuchten–Mualem soil hydraulic parameters.

This repository contains:
- A **Python script** (`pdt_function.py`) that evaluates the polynomial pedotransfer function.
- Precomputed **polynomial coefficients** (`Full_Polynomial.*`, `0.5_Polynomial.*`).
- A **web-based polynomial calculator** (GitHub Pages) for interactive use.


## Web calculator

You can use the polynomial calculator directly in the browser:
https://infoleon.github.io/Pedotransfer/

The Python script requires Sympy package to run. To install Sympy, please, activate the desired environment and run the command:

To install SymPy in your preferred environment:

```bash
pip install sympy


The script in this folder is a Python code to generate and solve the polynomial pedotransfer function associated to article https://doi.org/10.1016/j.geoderma.2021.115308
