from matplotlib import axes, pyplot as plt
from read_wfc import parse_file_cube
from read_wfc import parse_file_cube_atoms
import numpy as np

cube_URL = "./wfc/nvcenter1/spinup128.cube"

orbital_data, dd, nn = parse_file_cube(cube_URL)
atom_num, atoms_info = parse_file_cube_atoms(cube_URL)

rot_vec = np.array([1,1,1])
rot_vec = rot_vec / np.linalg.norm(rot_vec)
c3_theta = 2*np.pi/3
c3p_theta = 4*np.pi/3
carbon1 = np.array([2.67507, 0.88132, 0.88132])
carbon2 = np.array([0.88132, 2.67507, 0.88132])
carbon3 = np.array([0.88132, 0.88132, 2.67507])
ref_vec1 = 2*carbon1-carbon2-carbon3
ref_vec1 = ref_vec1 / np.linalg.norm(ref_vec1)
ref_vec2 = 2*carbon2-carbon1-carbon3
ref_vec2 = ref_vec2 / np.linalg.norm(ref_vec2)
ref_vec3 = 2*carbon3-carbon2-carbon3
ref_vec3 = ref_vec3 / np.linalg.norm(ref_vec3)

def rot_in_3d(k, v, theta):
    return v * np.cos(theta)+np.cross(k, v) * np.sin(theta)+k*np.dot(k, v)*(1-np.cos(theta))

def ref_in_3d(r, v):
    v_para = np.dot(r, v) * r
    v_vert = v - v_para
    return v_vert - v_para

def rotref_in_3d(r, k, v, theta):
    v = v * np.cos(theta)+np.cross(k, v) * np.sin(theta)+k*np.dot(k, v)*(1-np.cos(theta))
    v_para = np.dot(r, v) * r
    v_vert = v - v_para
    return v_vert - v_para

# def rotref_in_3d(r, k, v, theta):
#     return v * np.cos(theta)+np.cross(k, v) * np.sin(theta)+k*np.dot(k, v)*(1-np.cos(theta))

# v = np.array([2.67507,   0.88132,   0.88132])
# v = np.array([1.665450, 1.665450, 5.055144])
# vp = rot_in_3d(rot_vec, v, c3_theta)
# print(vp)
# print(dd)
# print(atoms_info)
# print(ref_in_3d(np.array([1, 0, 0]), np.array([1, 0, 0])))
# print(ref_vec1)
# print(ref_in_3d(ref_vec1, ref_vec1))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.quiver(0, 0, 0, ref_vec1[0], ref_vec1[1], ref_vec1[2], color='r', arrow_length_ratio=0.1)
# ax.quiver(0, 0, 0, v[0], v[1], v[2], color='b', arrow_length_ratio=0.1)
# ax.quiver(0, 0, 0, ref_in_3d(ref_vec1, v)[0], ref_in_3d(ref_vec1, v)[1], ref_in_3d(ref_vec1, v)[2], color='g', arrow_length_ratio=0.1)
# plt.show()

#=============================================

diag = {
    'e': [],
    'c3': [],
    'c3p': [],
    'sigma': [],
    'sigmap': [],
    'sigmapp': []
}
threshold = 1e-6
print(np.array(orbital_data).shape)
for ix, x in enumerate(orbital_data):
    for iy, y in enumerate(x):
        for iz, z in enumerate(y):
            if True:
                v = np.array([ix*dd[0], iy*dd[1], iz*dd[2]])
                vp_e = v
                vp_c3 = rot_in_3d(rot_vec, v, c3_theta)
                vp_c3p = rot_in_3d(rot_vec, v, c3p_theta)
                vp_sigma = ref_in_3d(ref_vec1, v)
                vp_sigmap = rotref_in_3d(ref_vec1, rot_vec, v, c3_theta)
                vp_sigmapp = rotref_in_3d(ref_vec1, rot_vec, v, c3p_theta)
                diag['e'].append(z * np.dot(v, vp_e))
                diag['c3'].append(z * np.dot(v, vp_c3))
                diag['c3p'].append(z * np.dot(v, vp_c3p))
                diag['sigma'].append(z * np.dot(v, vp_sigma))
                diag['sigmap'].append(z * np.dot(v, vp_sigmap))
                diag['sigmapp'].append(z * np.dot(v, vp_sigmapp))
                # if ix % 8 == 0 and iy % 8 == 0 and iz % 8 == 0:
                #     ax.quiver(v[0], v[1], v[2], vp_c3[0], vp_c3[1], vp_c3-[2], color='r', arrow_length_ratio=0.1)
            else:
                diag['e'].append(0)
                diag['c3'].append(0)
                diag['c3p'].append(0)
                diag['sigma'].append(0)
                diag['sigmap'].append(0)
                diag['sigmapp'].append(0)
            # print(ix*dd[0], iy*dd[1], iz*dd[2], z)
            
# plt.show()

red_ch_e = np.sum(np.array(diag['e']))
red_ch_c3 = np.sum(np.array(diag['c3']))
red_ch_c3p = np.sum(np.array(diag['c3p']))
red_ch_sigma = np.sum(np.array(diag['sigma']))
red_ch_sigmap = np.sum(np.array(diag['sigmap']))
red_ch_sigmapp = np.sum(np.array(diag['sigmapp']))

# red_ch_e = 9.002417842769734
# red_ch_c3 = 6.817714879474221
# red_ch_c3p = 6.81771487947422
# red_ch_sigma = 7.519884888101658
# red_ch_sigmap = 7.618647130081336
# red_ch_sigmapp = 7.572013513043791

norm = red_ch_e
red_ch_e = red_ch_e / norm
red_ch_c3 = red_ch_c3 / norm
red_ch_c3p = red_ch_c3p / norm
red_ch_sigma = red_ch_sigma / norm
red_ch_sigmap = red_ch_sigmap / norm
red_ch_sigmapp = red_ch_sigmapp / norm

print(red_ch_e)
print(red_ch_c3)
print(red_ch_c3p)
print(red_ch_sigma)
print(red_ch_sigmap)
print(red_ch_sigmapp)

a1 = red_ch_e*1+red_ch_c3*1+red_ch_c3p*1+red_ch_sigma*1+red_ch_sigmap*1+red_ch_sigmapp*1
a2 = red_ch_e*1+red_ch_c3*1+red_ch_c3p*1+red_ch_sigma*(-1)+red_ch_sigmap*(-1)+red_ch_sigmapp*(-1)
e = red_ch_e*2+red_ch_c3*1+red_ch_c3p*1+red_ch_sigma*(0)+red_ch_sigmap*(0)+red_ch_sigmapp*(0)

print("a1: ", a1)
print("a2: ", a2)
print("e: ", e)