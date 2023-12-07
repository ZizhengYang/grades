import numpy as np
from scipy.linalg import block_diag
from scipy.ndimage import affine_transform
from ase.io.cube import read_cube_data
from ase.io.cube import write_cube
from ase.io.xsf import read_xsf
from scipy.ndimage import affine_transform
from scipy.spatial.transform import Rotation
from scipy.linalg import fractional_matrix_power
# from westpy import VData

import math
import point_group
import read_wfc
from read_wfc import parse_file_cube

from symm import Cell
import sys

# cube_URL = "./wfc/nvcenter1/spindown124.cube"


def run_symmetry_analysis(URL, symmetry_group):
    VC3, VC2_12, VC2_23, VC2_31, Nor_12, Nor_23, Nor_31, origin, cell = point_group.find_cell_info(cube_URL)
    
    # cube_URL = "./wfc/nvcenter1/spindown124.cube"
    # cube_URL_2 = '/Users/zzyang/Documents/Pysymmetry/wfcs/nvcenter/spindown124.cube'

    #First read in any .cube data:
    # Testdata, atoms = read_cube_data('/Users/zzyang/Documents/Pysymmetry/wfcs/hbn/Banddn137.cube')
    # f = open('/Users/zzyang/Documents/Pysymmetry/wfcs/nvcenter/spindown124.xsf', "r")
    Testdata, atoms = read_cube_data(cube_URL)
    ase_atoms = atoms

    pg_name="C3v",
    # pg_operations = point_group.C3v(VC3, Nor_12, Nor_23, Nor_31, origin, cell)
    # pg_ctable = point_group.C3v_character(VC3, Nor_12, Nor_23, Nor_31, origin, cell)
    if degenerate == "False":
        pg_operations = point_group.PointGroup[symmetry_group](VC3, VC2_12, VC2_23, VC2_31, Nor_12, Nor_23, Nor_31, origin, cell)
        pg_ctable = point_group.PointGroupCharacter[symmetry_group]()
    elif degenerate == "True":
        pg_operations = point_group.PointGroup[symmetry_group](VC3, VC2_12, VC2_23, VC2_31, Nor_12, Nor_23, Nor_31, origin, cell)
        pg_ctable = point_group.PointGroupCharacter_degenerate[symmetry_group]()
    else:
        print("There is a problem with input .... ....")
    # print(pg_ctable)
    # pg_operations=point_group.C3v
    # pg_ctable=point_group.C3v_character

    #------------------------------------------------------------------------------------------------------------------------

    orbital_data, dd, nn = parse_file_cube(cube_URL)
    # orbital_data = VData(cube_URL, normalize="sqrt").data
    orbital_data = np.array(orbital_data)
    orbital_data = np.sign(orbital_data) * np.sqrt(np.abs(orbital_data))
    n_orb = 1
    omega = Cell(ase_atoms).omega
    # vdata.cell.omega

    N = nn[0] * nn[1] * nn[2]

    number_of_particle = 1
    assert number_of_particle >= 1

    rep_matrices_numlt = {}
    # for R1 in pg_operations:
    #     for R2 in pg_operations:
    #         rep_matrices.update({R1+" "+R2: np.zeros([n_orb, n_orb])})
    rep_matrices = {
        R: np.zeros([n_orb, n_orb]) for R in pg_operations
    }
    # print(orbital_data)
    for R, op in pg_operations.items():
        for j in range(n_orb):
            Rfj = op(orbital_data)
            for i in range(n_orb):
                fi = orbital_data
                print("Rfj is", np.array(Rfj).shape)
                print("fi is", np.array(fi).shape)
                rep_matrices[R][i, j] = omega / N * np.sum(fi * Rfj)
        print(rep_matrices[R])



if __name__ == "__main__":
    num_argv = len(sys.argv)
    symmetry = sys.argv[1]
    degenerate = sys.argv[2]
    for i in range(3, num_argv):
        cube_URL = sys.argv[i]
        run_symmetry_analysis(cube_URL, symmetry)
        
        
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/method3.py C3v True ./wfc/nvcenter1/spindown124.cube ./wfc/nvcenter1/spindown125.cube ./wfc/nvcenter1/spindown126.cube ./wfc/nvcenter1/spindown127.cube ./wfc/nvcenter1/spindown128.cube ./wfc/nvcenter1/spindown129.cube ./wfc/nvcenter1/spindown130.cube

# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py C3v False ./wfc/nvcenter1/spindown124.cube ./wfc/nvcenter1/spindown125.cube ./wfc/nvcenter1/spindown126.cube ./wfc/nvcenter1/spindown127.cube ./wfc/nvcenter1/spindown128.cube ./wfc/nvcenter1/spindown129.cube ./wfc/nvcenter1/spindown130.cube

# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Banddn138.cube wfc/hbn/Banddn139.cube wfc/hbn/Banddn140.cube wfc/hbn/Banddn141.cube wfc/hbn/Banddn142.cube wfc/hbn/Banddn143.cube wfc/hbn/Banddn144.cube wfc/hbn/Banddn145.cube wfc/hbn/Banddn146.cube wfc/hbn/Banddn147.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Bandup138.cube wfc/hbn/Bandup139.cube wfc/hbn/Bandup140.cube wfc/hbn/Bandup141.cube wfc/hbn/Bandup142.cube wfc/hbn/Bandup143.cube wfc/hbn/Bandup144.cube wfc/hbn/Bandup145.cube wfc/hbn/Bandup146.cube wfc/hbn/Bandup147.cube

# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h False wfc/hbn/Banddn139.cube wfc/hbn/Banddn140.cube wfc/hbn/Banddn141.cube wfc/hbn/Banddn142.cube wfc/hbn/Banddn143.cube wfc/hbn/Banddn144.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Banddn139.cube wfc/hbn/Banddn140.cube wfc/hbn/Banddn141.cube wfc/hbn/Banddn142.cube wfc/hbn/Banddn143.cube wfc/hbn/Banddn144.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h False wfc/hbn/Bandup139.cube wfc/hbn/Bandup140.cube wfc/hbn/Bandup141.cube wfc/hbn/Bandup142.cube wfc/hbn/Bandup143.cube wfc/hbn/Bandup144.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Bandup139.cube wfc/hbn/Bandup140.cube wfc/hbn/Bandup141.cube wfc/hbn/Bandup142.cube wfc/hbn/Bandup143.cube wfc/hbn/Bandup144.cube