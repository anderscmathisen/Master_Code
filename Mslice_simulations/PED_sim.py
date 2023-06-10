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


def uniform_circle(N, radius=5):
    thethas = np.linspace(0, 360, N + 1)

    xs = radius * np.cos(np.deg2rad(thethas))
    ys = radius * np.sin(np.deg2rad(thethas))

    arr = np.array([xs, ys]).T

    return arr[:-1]


def shift_image(X, dx, dy):
    X = np.roll(X, dy, axis=1)
    X = np.roll(X, dx, axis=2)
    if dy > 0:
        X[:, :dy, :] = 0
    elif dy < 0:
        X[:, dy:, :] = 0
    if dx > 0:
        X[:, :, :dx] = 0
    elif dx < 0:
        X[:, :, dx:] = 0
    return X


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


Num_beam_tilts = 100
circ_radius = 17.453293  # 1 deg precession

out_arr = np.zeros((300, 1365, 1365))

#loop over the N beam tilts
for tilt in tqdm(uniform_circle(Num_beam_tilts, radius=circ_radius)):

    output = pyms.CBED(
        crystal,
        gridshape,
        eV,
        app,
        thicknesses,
        beam_tilt=tilt,
        tilt_units="mrad",
        subslices=subslices,
        tiling=tiling,
        nfph=nfph,
        showProgress=False,
        device_type=device,
    )

    output = np.array(output)

    k = wavev(eV)
    tilt_ = np.asarray(tilt) * 1e-3 * k

    shift = tilt_ * crystal.unitcell[:2] * np.asarray(tiling)

    x, y = shift
    #Shift imagestack to center zero beam
    output = shift_image(output, int(-y), int(-x))
    out_arr += output


#dump results with pickle
with open("output_full.pickle", "wb") as f:
    pickle.dump(out_arr, f)
