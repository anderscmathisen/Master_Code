import pyms
import numpy as np
import torch
import matplotlib.pyplot as plt
from pyms.Probe import wavev
from tqdm import tqdm
import pickle


crystal = pyms.structure.fromfile(
    "./data/Rhom_ErMnO3_Ppos_cba.p1",
    temperature_factor_units="B",
    atomic_coordinates="fractional",
)



# Subslicing of crystal for multislice
slices = 2
subslices = [(i + 1) / slices for i in range(slices)]


# Grid size in pixels
gridshape = [2048, 2048]

# Tile out grid to acceptable size
tiling = [32, 32]

# Probe accelerating voltage in eV
eV = 200e3

# Probe forming aperture in mrad
app = 1.23121  # gives about same radius as experimental

# Number of frozen phonon passes
nfph = 25  # the default value


device = torch.device("cpu")

thicknesses = np.linspace(250, 1500, 300)

out_arr = np.zeros((300, 1365, 1365))


output = pyms.CBED(
        crystal,
        gridshape,
        eV,
        app,
        thicknesses,
        subslices=subslices,
        tiling=tiling,
        nfph=nfph,
        showProgress=False,
        device_type=device,
    )


#dump results with pickle
with open("output_full_CBED.pickle", "wb") as f:
    pickle.dump(out_arr, f)
