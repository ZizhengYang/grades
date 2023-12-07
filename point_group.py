from os import urandom
import numpy as np
from scipy.linalg import block_diag
from scipy.ndimage import affine_transform
from ase.io.cube import read_cube_data
from ase.io.cube import write_cube
from ase.io.xsf import read_xsf
from scipy.ndimage import affine_transform
from scipy.spatial.transform import Rotation

import math
import read_wfc
from read_wfc import parse_file_cube

from symm import Cell

class PointGroupOp:
    
    def __init__(self, T:np.ndarray, origin=None, cell:np.ndarray=None):
        assert T.shape == (4, 4)
        self.T = T
        if cell is not None:
            assert np.shape(cell)==(3,3)
            if origin is not None:
                origin = np.linalg.inv(cell) @ origin
            self.set_coord(cell)
        if origin is not None:
            assert len(origin) == 3
            self.set_origin(origin)
            
    def set_coord(self, cell):
        TC = block_diag(cell, 1)
        self.T = np.linalg.inv(TC) @ self.T @ TC
        
    def set_origin(self, origin):
        x0, y0, z0 = origin
        TR = np.array(
            [
                [1, 0, 0, x0],
                [0, 1, 0, y0],
                [0, 0, 1, z0],
                [0, 0, 0, 1],
            ]
        )
        self.T = TR @ self.T @ np.linalg.inv(TR)
        
    @property
    def inv(self):
        return PointGroupOp(T=np.linalg.inv(self.T))
    
    def __matmul__(self, other):
        assert isinstance(other, PointGroupOp)
        return PointGroupOp(T=self.T @ other.T)
    
    def __call__(self, f):
        return affine_transform(f, matrix=self.T, mode='grid-wrap')


class PointGroupRef(PointGroupOp):
    def __init__(
        self,
        normal,
        origin = (0.0, 0.0, 0.0),
        cell: np.ndarray = None
    ):
        a, b, c = np.array(normal) / np.linalg.norm(normal)
        RE = np.array(
            [
                [1 - 2 * a * a, -2 * a * b, -2 * a * c, 0],
                [-2 * a * b, 1 - 2 * b * b, -2 * b * c, 0],
                [-2 * a * c, -2 * b * c, 1 - 2 * c * c, 0],
                [0, 0, 0, 1],
            ]
        )
        super(PointGroupRef, self).__init__(T=RE, origin=origin,cell=cell)


class PointGroupRot(PointGroupOp):
    def __init__(
        self,
        rotvec,
        origin = (0.0, 0.0, 0.0),
        cell: np.ndarray = None
    ):
        rotation = Rotation.from_rotvec(rotvec)
        RO = block_diag(rotation.as_matrix().T, 1)
        super(PointGroupRot, self).__init__(T=RO, origin=origin,cell=cell)


class PointGroupInv(PointGroupOp):
    def __init__(self, origin = (0.0, 0.0, 0.0),
                 cell:np.ndarray = None):
        super(PointGroupInv, self).__init__(
            T=block_diag(-1 * np.eye(3), 1), origin=origin, cell=cell
        )

class PointGroupRotRef(PointGroupOp):
    def __init__(
        self,
        rotvec,
        origin = (0.0, 0.0, 0.0),
        multiple: int = 1 ,
        cell: np.ndarray = None
    ):
        R0 = Rotation.from_rotvec(rotvec).as_matrix()
        a, b, c = np.array(rotvec) / np.linalg.norm(rotvec)
        RE = np.array(
            [
                [1 - 2 * a * a, -2 * a * b, -2 * a * c],
                [-2 * a * b, 1 - 2 * b * b, -2 * b * c],
                [-2 * a * c, -2 * b * c, 1 - 2 * c * c],
            ]
        )
        S = RE @ R0
        if multiple > 1:
            for m in range(1, multiple):
                S = RE@R0@S
        S_Affine = block_diag(S.T, 1)
        super(PointGroupRotRef, self).__init__(T=S_Affine, origin=origin, cell=cell)


def find_cell_info(cube_URL):# -> tuple[NDArray[Any], Any, Any, Any, Any, Any]:
    #First read in any .cube data:
    # Testdata, atoms = read_cube_data('/Users/zzyang/Documents/Pysymmetry/wfcs/hbn/Banddn137.cube')
    # f = open('/Users/zzyang/Documents/Pysymmetry/wfcs/nvcenter/spindown124.xsf', "r")
    Testdata, atoms = read_cube_data(cube_URL)
    ase_atoms = atoms
    read_wfc.parse_file_cube_atoms(cube_URL)
    # Testdata = parse_file_xsf('/Users/zzyang/Documents/Pysymmetry/wfcs/nvcenter/spindown124.xsf')
    # print(Testdata)
    # print(atoms)
    # print("Prefactor of Origin determined by real-space grid:",np.shape(Testdata))
    #origin = Origin*np.array(np.shape(Testwfc.data)) #origin is in the Vn coordinate
    #The Cartesian coordinate of three N atoms around the vacancy
    A=atoms.cell[:].T #This enable A x_new=n_old transform
    # print(atoms.cell)
    # print(A)
    # N1=A@np.array([0.38594,   0.61406,   0.00000])
    # N2=A@np.array([0.38594,   0.43855,   0.00000])
    # N3=A@np.array([0.56145,   0.61406,   0.00000])
    N1=A@np.array([2.67507,   0.88132,   0.88132])
    N2=A@np.array([0.88132,   2.67507,   0.88132])
    N3=A@np.array([0.88132,   0.88132,   2.67507])
    #Origin point is :
    origin = (N1+N2+N3)/3
    #C3,S3 Rotation axis(normalized); also the sigma_h's normal vector for rotation plane
    VC3=np.array([0,0,1])
    #C2 rotation axis(normalized) :
    VC2_12 = (N3-N1)+(N3-N2) ; VC2_12 = VC2_12/np.linalg.norm(VC2_12)
    VC2_23 = (N1-N3)+(N1-N2) ; VC2_23 = VC2_23/np.linalg.norm(VC2_23)
    VC2_31 = (N2-N1)+(N3-N1) ; VC2_31 = VC2_31/np.linalg.norm(VC2_31)
    #Sigma_v rotation plane 's Normal vector
    Nor_12 = np.cross(VC3, VC2_12)
    Nor_23 = np.cross(VC3, VC2_23)
    Nor_31 = np.cross(VC3, VC2_31)
    [a,b,c] = 1/np.array(np.shape(Testdata))
    T = block_diag(a,b,c) #This is the scale factor
    cell = A@T #M@v_n=v_xy
    return VC3, VC2_12, VC2_23, VC2_31, Nor_12, Nor_23, Nor_31, origin, cell



def C3v(VC3, VC2_12, VC2_23, VC2_31, Nor_12, Nor_23, Nor_31, origin, cell):
    return {
            "E": PointGroupOp(T=np.eye(4)),
            "C3_1": PointGroupRot(rotvec=2 * np.pi / 3 * VC3, origin=origin, cell = cell),
            "C3_2": PointGroupRot(rotvec=4 * np.pi / 3 * VC3, origin=origin, cell = cell),
            "Sig_v1" : PointGroupRef(normal=Nor_12,origin=origin, cell = cell),
            "Sig_v2" : PointGroupRef(normal=Nor_23,origin=origin, cell = cell),
            "Sig_v3" : PointGroupRef(normal=Nor_31,origin=origin, cell = cell)
    }
    
def D3h(VC3, VC2_12, VC2_23, VC2_31, Nor_12, Nor_23, Nor_31, origin, cell):
    return {
        "E": PointGroupOp(T=np.eye(4)),
        "C3_1": PointGroupRot(rotvec=2 * np.pi / 3 * VC3,
                                   origin=origin, cell = cell),
        "C3_2": PointGroupRot(rotvec=4 * np.pi / 3 * VC3,
                                   origin=origin, cell = cell),
        "C2_1" : PointGroupRot(rotvec=np.pi*VC2_12,
                                   origin=origin, cell = cell),
        "C2_2" : PointGroupRot(rotvec=np.pi*VC2_23,
                                   origin=origin, cell = cell),
        "C2_3" : PointGroupRot(rotvec=np.pi*VC2_31,
                                   origin=origin, cell = cell),
        "Sig_h": PointGroupRef(normal=VC3, origin=origin, cell = cell),
        "S3_1" : PointGroupRotRef(rotvec=2 * np.pi / 3 * VC3,
                                           origin=origin, multiple=1, cell = cell),
        "S3_2" : PointGroupRotRef(rotvec=2 * np.pi / 3 * VC3,
                                           origin=origin, multiple=2, cell = cell),
        "Sig_v1" : PointGroupRef(normal=Nor_12,origin=origin, cell = cell),
        "Sig_v2" : PointGroupRef(normal=Nor_23,origin=origin, cell = cell),
        "Sig_v3" : PointGroupRef(normal=Nor_31,origin=origin, cell = cell)
    }

def C3v_character():
    return {
            "A1": [1,1,1,1,1,1],
            "A2": [1,1,1,-1,-1,-1],
            "E": [2,-1,-1,0,0,0]
    }

def C3v_character_degenerate():
    return {
            "A1": [1,1,1,1,1,1],
            "A2": [1,1,1,-1,-1,-1],
            "Ex": [1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2,1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2],
            "Ey": [1,-1/2+math.sqrt(3)/2,-1/2-math.sqrt(3)/2,-1,1/2-math.sqrt(3)/2,1/2+math.sqrt(3)/2]
    }

def D3h_character():
    return {
        "A1p": [1,1,1,1,1,1,1,1,1,1,1,1],
        "A2p": [1,1,1,-1,-1,-1,1,1,1,-1,-1,-1],
        "Ep": [2,-1,-1,0,0,0,2,-1,-1,0,0,0],
        "A1pp" : [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
        "A2pp" : [1,1,1,-1,-1,-1,-1,-1,-1,1,1,1],
        "Epp" : [2,-1,-1,0,0,0,-2,1,1,0,0,0]
    }
    

def D3h_character_degenerate():
    return {
        "A1p": [1,1,1,1,1,1,1,1,1,1,1,1],
        "A2p": [1,1,1,-1,-1,-1,1,1,1,-1,-1,-1],
        "Epx": [1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2,1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2,1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2,1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2],
        "Epy": [1,-1/2+math.sqrt(3)/2,-1/2-math.sqrt(3)/2,-1,1/2-math.sqrt(3)/2,1/2+math.sqrt(3)/2,1,-1/2+math.sqrt(3)/2,-1/2-math.sqrt(3)/2,-1,1/2-math.sqrt(3)/2,1/2+math.sqrt(3)/2],
        "A1pp" : [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
        "A2pp" : [1,1,1,-1,-1,-1,-1,-1,-1,1,1,1],
        "Eppx" : [1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2,1,-1/2-math.sqrt(3)/2,-1/2+math.sqrt(3)/2,-1,-(-1/2-math.sqrt(3)/2),-(-1/2+math.sqrt(3)/2),-1,-(-1/2-math.sqrt(3)/2),-(-1/2+math.sqrt(3)/2)],
        "Eppy" : [1,-1/2+math.sqrt(3)/2,-1/2-math.sqrt(3)/2,-1,1/2-math.sqrt(3)/2,1/2+math.sqrt(3)/2,-1,-(-1/2+math.sqrt(3)/2),-(-1/2-math.sqrt(3)/2),1,-(1/2-math.sqrt(3)/2),-(1/2+math.sqrt(3)/2)],
    }
    

PointGroup = {
    "C3v": C3v,
    "D3h": D3h,
}

PointGroupCharacter = {
    "C3v": C3v_character,
    "D3h": D3h_character,
}

PointGroupCharacter_degenerate = {
    "C3v": C3v_character_degenerate,
    "D3h": D3h_character_degenerate,
}