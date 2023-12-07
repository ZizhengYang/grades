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
        # print(rep_matrices[R])
    # for R1, op1 in rep_matrices.items():
    #     for R2, op2 in rep_matrices.items():
    #         rep_matrices_numlt.update({R1+" "+R2: op1 @ op2})
    # for R, D in rep_matrices.items():
    #     print(R)
    for R, D in rep_matrices.items():
        S = D @ D.T
        U = fractional_matrix_power(S, -1 / 2)
        D[...] = U @ D
    # print("+"*15)
    # for k, v in rep_matrices.items():
    #     print(k)
    #     print(v)
    # print("+"*15)
    # for k, v in rep_matrices_numlt.items():
    #     print(k)
    #     print(v)
    # print("+"*15)
        
    symms = []

    # for R in rep_matrices.values():
    #     print(R)
    vec = np.zeros(n_orb)
    vec[0] = 1
    irprojs = {}
    h = len(pg_operations)
    for irrep, chis in pg_ctable.items():
        l = chis[0]
        pvec = np.zeros_like(vec)
        for chi, U in zip(chis, rep_matrices.values()):
            pvec += chi * U @ vec
        irprojs[irrep] = l / h * np.sum(vec * pvec)
        print("vec is: ", vec)
        print("pvec is: ", pvec)
    irreps = list(irprojs.keys())
    irproj_values = list(irprojs.values())
    # print(irproj_values)
    imax = np.argmax(irproj_values)
    # print(irprojs)
    symms.append(f"{irreps[imax]}({irproj_values[imax]:.6f})")

    print("="*50)
    print("Apply symmetry analysis to:", URL)
    print("The final results are:", irprojs)
    print("Irrep of orbitals:", symms)

# multi_particle_irprojs = {}
# for R1, p1 in irprojs.items():
#     for R2, p2 in irprojs.items():
#         # print(p1)
#         multi_particle_irprojs.update({R1+" "+R2: p1 * p2})
#         print(R1+" "+R2)
#         print(p1 * p2)
# imax = np.argmax(multi_particle_irprojs)
# multi_particle_irprojs_values = list(multi_particle_irprojs.values())
# multi_particle_irprojs_keys = list(multi_particle_irprojs.keys())
# symms.append(f"{multi_particle_irprojs_keys[imax]}({multi_particle_irprojs_values[imax]:.2f})")
# print("Irrep of orbitals:", symms)


if __name__ == "__main__":
    num_argv = len(sys.argv)
    symmetry = sys.argv[1]
    degenerate = sys.argv[2]
    for i in range(3, num_argv):
        cube_URL = sys.argv[i]
        run_symmetry_analysis(cube_URL, symmetry)
        
        
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py C3v True ./wfc/nvcenter1/spindown124.cube ./wfc/nvcenter1/spindown125.cube ./wfc/nvcenter1/spindown126.cube ./wfc/nvcenter1/spindown127.cube ./wfc/nvcenter1/spindown128.cube ./wfc/nvcenter1/spindown129.cube ./wfc/nvcenter1/spindown130.cube

# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py C3v False ./wfc/nvcenter1/spindown124.cube ./wfc/nvcenter1/spindown125.cube ./wfc/nvcenter1/spindown126.cube ./wfc/nvcenter1/spindown127.cube ./wfc/nvcenter1/spindown128.cube ./wfc/nvcenter1/spindown129.cube ./wfc/nvcenter1/spindown130.cube

# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Banddn138.cube wfc/hbn/Banddn139.cube wfc/hbn/Banddn140.cube wfc/hbn/Banddn141.cube wfc/hbn/Banddn142.cube wfc/hbn/Banddn143.cube wfc/hbn/Banddn144.cube wfc/hbn/Banddn145.cube wfc/hbn/Banddn146.cube wfc/hbn/Banddn147.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Bandup138.cube wfc/hbn/Bandup139.cube wfc/hbn/Bandup140.cube wfc/hbn/Bandup141.cube wfc/hbn/Bandup142.cube wfc/hbn/Bandup143.cube wfc/hbn/Bandup144.cube wfc/hbn/Bandup145.cube wfc/hbn/Bandup146.cube wfc/hbn/Bandup147.cube

# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h False wfc/hbn/Banddn139.cube wfc/hbn/Banddn140.cube wfc/hbn/Banddn141.cube wfc/hbn/Banddn142.cube wfc/hbn/Banddn143.cube wfc/hbn/Banddn144.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Banddn139.cube wfc/hbn/Banddn140.cube wfc/hbn/Banddn141.cube wfc/hbn/Banddn142.cube wfc/hbn/Banddn143.cube wfc/hbn/Banddn144.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h False wfc/hbn/Bandup139.cube wfc/hbn/Bandup140.cube wfc/hbn/Bandup141.cube wfc/hbn/Bandup142.cube wfc/hbn/Bandup143.cube wfc/hbn/Bandup144.cube
# /Users/zzyang/opt/anaconda3/envs/my_research/bin/python3.11 /Users/zzyang/Documents/grades/grades/main.py D3h True wfc/hbn/Bandup139.cube wfc/hbn/Bandup140.cube wfc/hbn/Bandup141.cube wfc/hbn/Bandup142.cube wfc/hbn/Bandup143.cube wfc/hbn/Bandup144.cube