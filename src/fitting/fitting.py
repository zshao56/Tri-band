import matplotlib.pylab as plt
import numpy as np
import tidy3d as td
from tidy3d.plugins.dispersion import FastDispersionFitter, AdvancedFastFitterParam

# List of materials
# materials = ['ITO-BaF2', 'ITO-glass', 'SPE', "WO3-coloration", "WO3-bleach"]  # Add more materials as needed
materials = ['nk/WO3-bleach-modified'] 
# Loop over each material
for material in materials:
    print(f"Processing {material}...")

    # File names
    fname = material + ".csv"
    fname1 = material + ".json"

    # Load nk data from CSV file
    fitter = FastDispersionFitter.from_file(fname, skiprows=1, delimiter=",", encoding="utf-8")



    # Advanced fitting parameters
    advanced_param = AdvancedFastFitterParam(weights=(1, 0.5))

    # Fit the model
    medium, rms_error = fitter.fit(max_num_poles=8, advanced_param=advanced_param, tolerance_rms=2e-2)

    # Update wavelength range and plot fitted data
    fitter = fitter.copy(update={"wvl_range": (0.4, 20)})
    fitter.plot(medium)
    plt.title(f"{material} - Fitted Data")
    plt.savefig(f'{material}_fitting.jpg')

    # Save poles to JSON file
    medium.to_file(fname1)
    print(f"Saved poles to {fname1}")

    # Optional: Load the file back to verify
    medium = td.PoleResidue.from_file(fname1)
    print(f"Loaded medium from {fname1}\n")

print("Processing complete for all materials.")
