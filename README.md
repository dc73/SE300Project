# SE300 FEA Project

A desktop demonstration that converts an STL surface mesh to the project's
text input format and visualizes placeholder result fields.

> The current `run_analysis` function does **not** solve finite-element
> equations. Its values are deterministic placeholders for testing the GUI.

## Setup and run

Requires Python 3.10+ and Tk.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python appFEA.py
```

## Layout

- `fea/` — GUI, conversion-file I/O, materials, and analysis reader
- `examples/` — standalone experiments and the original MATLAB comparison
- `docs/` — project demo media
- `tests/` — fast checks for the input/output pipeline

Run checks with `python -m unittest`.
