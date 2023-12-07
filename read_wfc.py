import re
import numpy as np
# import ase.io.cube
# from westpy import VData
# from westpy import Angstrom
import const

def parse_file_cube_atoms(path):
    f = open(path, "r")
    lines = f.readlines()
    line_atom = lines[2]
    atom_num = int(remove_space_keep_num(line_atom)[0])
    atoms_info = []
    for i in range(atom_num):
        atom_info = remove_space_keep_num(lines[i+6])
        # print(remove_space_keep_num(lines[i+6]))
        atom_info[0] = float(atom_info[0])
        atom_info[1] = const.atom_number[atom_info[0]]
        for j in range(2, 5):
            atom_info[j] = float(atom_info[j])
        # print(atom_info)
        atoms_info.append(atom_info)
    return atom_num, atoms_info


def parse_file_cube(path):
    f = open(path, "r")
    lines = f.readlines()
    line_atom = lines[2]
    atom_num = int(remove_space_keep_num(line_atom)[0])
    line_x = lines[3]
    line_y = lines[4]
    line_z = lines[5]
    param_x = remove_space_keep_num(line_x)
    param_y = remove_space_keep_num(line_y)
    param_z = remove_space_keep_num(line_z)
    nx = int(param_x[0])
    ny = int(param_y[0])
    nz = int(param_z[0])
    dx = float(param_x[1])
    dy = float(param_y[2])
    dz = float(param_z[3])
    data = []
    for i in range(atom_num+6, len(lines)):
        data_i = remove_space_keep_num(lines[i])
        for data_ij in data_i:
            data.append(float(data_ij))
    data_x = []
    count = 0
    for x in range(nx):
        data_y = []
        for y in range(ny):
            data_z = []
            for z in range(nz):
                data_z.append(data[count])
                count += 1
            data_y.append(data_z)
        data_x.append(data_y)
    angstrom = const.angstrom
    return data_x, [dx/angstrom, dy/angstrom, dz/angstrom], [nx, ny, nz]


def remove_space_keep_num(in_str):
    str_params = re.split(r"[;,\s]\s*", in_str)
    ret_array = []
    for str_param in str_params:
        if not str_param == '':
            ret_array.append(str_param)
    return ret_array


# my_data, dd, nn = parse_file_cube("./wfc/nvcenter/spindown124.cube")
# west_data, _ = ase.io.cube.read_cube_data("./wfc/nvcenter/spindown124.cube")
# print(np.array(my_data).shape)
# print(np.array(west_data).shape)
# vdata = VData("./wfc/nvcenter/spindown124.cube", normalize='sqrt')
# print(np.array(my_data[0][0]))
# print(np.array(west_data[0][0]))
# print((my_data == west_data)[0][0])
# print(vdata.data == west_data)
# print(Angstrom)
# print(angstrom)
# print(vdata.dx)
# print(vdata.cell.R[(0,0)])
# print(dd)
# print(nn)

# atom_num, atoms_info = parse_file_cube_atoms("./wfc/nvcenter/spindown124.cube")
# print(atom_num)
# print(atoms_info)