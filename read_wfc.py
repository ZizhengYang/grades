import re
import numpy as np

from grades.const import dict_atom_number

class ParseCube:
    def __init__(self, file_path):
        """
        Initialize the ParseCube object with the path to the Gaussian Cube file.

        Parameters:
        - file_path (str): Path to the Gaussian Cube file.
        """
        self.file_path = file_path
        self.lines = []
        self.num_atoms = 0
        self.sampling_matrix = None
        self.sampling_step = []
        self.origin = (0.0, 0.0, 0.0)
        self.atom_coordinates = {}
        self.volume_data = []
        
    def remove_space_keep_num(self, in_str):
        str_params = re.split(r"[;,\s]\s*", in_str)
        ret_array = []
        for str_param in str_params:
            if not str_param == '':
                ret_array.append(str_param)
        return ret_array

    def read_cube_file(self):
        """
        Read the Gaussian Cube file and extract relevant information.

        This function reads the header information, atomic coordinates, and volume data.
        """
        with open(self.file_path, 'r') as cube_file:
            self.lines = cube_file.readlines()
            self.read_atom_numbers(cube_file)
            
            # Read and parse header information
            self.read_header_sampling_steps(cube_file)

            # Read atomic coordinates
            self.read_atom_coordinates(cube_file)

            # Read volume data
            # self.read_volume_data(cube_file)
            
    def read_atom_numbers(self, cube_file):
        """
        Reads the number of atoms from a Gaussian Cube file.

        Parameters:
        - cube_file (file): A file object representing the Gaussian Cube file.
        - self.num_atoms (int): Number of atoms in the molecular structure.
        """
        self.num_atoms = int(self.remove_space_keep_num(self.lines[2])[0])

    def read_header_sampling_steps(self, cube_file):
        """
        Read the header information from the Gaussian Cube file.

        Parameters:
        - cube_file (file): Open file object for the Gaussian Cube file.
        """
        # Read the number of atoms and volume dimensions
        param_x = self.remove_space_keep_num(self.lines[3])
        param_y = self.remove_space_keep_num(self.lines[4])
        param_z = self.remove_space_keep_num(self.lines[5])
        
        self.sampling_matrix = np.array([
            [float(param_x[1]), float(param_x[2]), float(param_x[3])],
            [float(param_y[1]), float(param_y[2]), float(param_y[3])],
            [float(param_z[1]), float(param_z[2]), float(param_z[3])]
        ])
        
        self.sampling_step.append(float(param_x[1]))
        self.sampling_step.append(float(param_y[2]))
        self.sampling_step.append(float(param_z[3]))

        # Read the origin
        # self.origin = tuple(map(float, cube_file.readline().split()[1:4]))

    def read_atom_coordinates(self, cube_file):
        """
        Read atomic coordinates from the Gaussian Cube file.

        Parameters:
        - cube_file (file): Open file object for the Gaussian Cube file.
        """
        # for _ in range(self.num_atoms):
        #     atom_info = tuple(map(float, cube_file.readline().split()[2:5]))
        #     self.atom_coordinates.append(atom_info)
        for i in range(6, self.num_atoms+6):
            element = dict_atom_number[int(self.remove_space_keep_num(self.lines[i])[0])]
            position = [float(pos) for pos in self.remove_space_keep_num(self.lines[i])[2:]]
            self.atom_coordinates[f'{element}{i-6}'] = position
            
    def find_system_origin(self):
        return 0

    def read_volume_data(self, cube_file):
        """
        Read volume data from the Gaussian Cube file.

        Parameters:
        - cube_file (file): Open file object for the Gaussian Cube file.
        """
        for i in range(self.num_atoms+6, len(self.lines)):
            data_i = self.remove_space_keep_num(self.lines[i])
            for data_ij in data_i:
                print(data_ij)

    def get_num_atoms(self):
        """
        Get the number of atoms in the Gaussian Cube file.

        Returns:
        - int: Number of atoms.
        """
        return self.num_atoms

    def get_sampling_step(self):
        """
        Get the dimensions of the volume grid in the Gaussian Cube file.

        Returns:
        - tuple: (nx, ny, nz), representing the number of grid points in each dimension.
        """
        return self.sampling_step

    def get_origin(self):
        """
        Get the origin coordinates of the volume grid in the Gaussian Cube file.

        Returns:
        - tuple: (ox, oy, oz), representing the origin coordinates.
        """
        return self.origin
    
    

cube_parser = ParseCube("./wfc/nvcenter/spindown124.cube")
cube_parser.read_cube_file()

num_atoms = cube_parser.get_num_atoms()
sampling_step = cube_parser.get_sampling_step()
origin = cube_parser.get_origin()

print("Number of Atoms:", num_atoms)
print("Volume Dimensions:", sampling_step)
# print("Atom Positions:", cube_parser.atom_coordinates)
print("Origin:", origin)