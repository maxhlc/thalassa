# Add path for PyTHALASSA
import os
import sys

sys.path.append(
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "lib",
        )
    )
)

# Import PyTHALASSA
import pythalassa

# Import Numpy
import numpy as np

# Declare model
model = pythalassa.Model()

# Declare paths
paths = pythalassa.Paths()
paths.phys_path = "./data/physical_constants.txt"
paths.earth_path = "./data/earth_potential/GRIM5-S1.txt"
paths.kernel_path = "./data/kernels_to_load.furnsh"
paths.eop_path = "./data/eop_data.txt"

# Declare settings
settings = pythalassa.Settings()

# Declare spacecraft
spacecraft = pythalassa.Spacecraft()
spacecraft.mass = 8500.0
spacecraft.area_drag = 13.0
spacecraft.area_srp = 13.0
spacecraft.cd = 2.2
spacecraft.cr = 1.5

# Declare initial state
times = np.linspace(59043.0, 59043.0 + 365.25, 367)
stateIn = np.array([7000.0, 0.0, 0.0, 0.0, 8.0, 0.0])

# Declare propagator
propagator = pythalassa.Propagator(model, paths, settings, spacecraft)

# Propagate with THALASSA
statesOut = propagator.propagate(times, stateIn)
