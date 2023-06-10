import numpy as np
import matplotlib.pyplot as plt
from diffsims.generators.rotation_list_generators import get_beam_directions_grid
import diffpy
from diffsims.libraries.structure_library import StructureLibrary
from diffsims.generators.diffraction_generator import DiffractionGenerator
from diffsims.generators.library_generator import DiffractionLibraryGenerator
from diffsims.libraries.diffraction_library import (
    DiffractionLibrary,
    load_DiffractionLibrary,
)
import os
from orix.quaternion import Orientation, Rotation, symmetry


def GetDiffLibrary(
    diffraction_calibration: float,
    camera_length: int,
    half_radius: int,
    resolution: int = 1,
    cif_file: str = "/home/anderscm/Documents/Data/2022-09-13T113412.38676Z-SED and SPED data 1 and .5 presecion angle-5435/ErMnO3.cif",
    grid_cub: str = None,
    make_new: bool = False,
    minimum_intensity: float = 0.2,
    max_excitation_error: float = 0.08,
    center_spot : bool = False,
    precession_angle = 1
) -> DiffractionLibrary:

    cif_material = os.path.basename(cif_file).split(".")[0]

    filename = f"{cif_material}-{half_radius}-{camera_length}-{resolution}.difflib"

    basepath = os.path.dirname(os.path.abspath(__file__))

    filepath = os.path.join(basepath, "Libraries", filename)

    if os.path.exists(filepath) and not make_new:
        # pass
        return load_DiffractionLibrary(filepath, True)

    reciprocal_radius = (
        np.sqrt(half_radius**2 + half_radius**2) * diffraction_calibration
    )

    # importing the structure: cif file
    structure_matrix = diffpy.structure.loadStructure(cif_file)

    # parameters that determine how the templates are calculated
    diff_gen = DiffractionGenerator(
        accelerating_voltage=200,
        precession_angle=precession_angle,
        scattering_params=None,#"lobato",
        shape_factor_model="linear",
        minimum_intensity=minimum_intensity,
    )

    lib_gen = DiffractionLibraryGenerator(diff_gen)
    if grid_cub is None:
        grid_cub = get_beam_directions_grid(
            "hexagonal", resolution, mesh="spherified_cube_edge"
        )

    # library containing structures and orientations
    library_phases = StructureLibrary(["ErMnO3"], [structure_matrix], [grid_cub])

    # final diffraction library
    diff_lib = lib_gen.get_diffraction_library(
        library_phases,
        calibration=diffraction_calibration,
        reciprocal_radius=reciprocal_radius,
        half_shape=half_radius,
        with_direct_beam=center_spot,
        max_excitation_error=max_excitation_error,
    )
    diff_lib.pickle_library(filepath)

    return diff_lib


if __name__ == "__main__":
    calib = 0.015920479420444417 / 2
    CL = 10
    Half_R = 64 * 2
    arr1 = []
    with open("rotmat.txt", "r") as f:
        for line in f.readlines():
            arr = line.split("  ")
            temp_arr = []
            for elm in arr:
                temp_arr.append(float(elm))
            arr1.append(temp_arr)
    matrix = np.array(arr1).T

    orientation = Orientation.from_matrix(matrix, symmetry.C6h)

    

    orientation.scatter(projection = "ipf")
    plt.show()
    euler = np.degrees(orientation.to_euler())
    grid = np.array([euler[0]])

    difflib = GetDiffLibrary(
        calib,
        CL,
        Half_R,
        resolution=5,
        grid_cub=grid,
        make_new=True,
        minimum_intensity=0.01,
        max_excitation_error=0.15,
        center_spot = True
    )["ErMnO3"]
    sims = difflib["simulations"]
    oris = difflib["orientations"]
    for sim in sims:
        sim.plot()
        plt.show()
        # plt.imshow(sim.get_diffraction_pattern(sigma = 5))
        # plt.show()
