import numpy as np

from grades.const import dict_atom_number
import util

class Structure():
    def __init__(self, file_path):
        self.num_atoms = 0
        self.sampling_matrix = None
        self.sampling_step = []
        self.origin = (0.0, 0.0, 0.0)
        self.atom_coordinates = {}
        self.volume_data = []


class Wavefunction():
    def __init__(self, file_path):
        self.file_path = file_path
        self._lines = []
        self.num_atoms = 0
        self.sampling_matrix = None
        self.sampling_step = []
        self.origin = (0.0, 0.0, 0.0)
        self.atom_coordinates = {}
        self.volume_data = []
        
    def parse_cube(self):
        with open(self.file_path, 'r') as cube_file:
            self._lines = cube_file.readlines()
            self.read_atom_numbers(cube_file)
            self.read_header_sampling_steps(cube_file)      # Read and parse header information
            self.read_atom_coordinates(cube_file)           # Read atomic coordinates
            
    def read_atom_numbers(self):
        self.num_atoms = int(util.remove_space_keep_num(self._lines[2])[0])
            
    def read_header_sampling_steps(self, cube_file):
        param_x = util.remove_space_keep_num(self._lines[3])
        param_y = util.remove_space_keep_num(self._lines[4])
        param_z = util.remove_space_keep_num(self._lines[5])
        
        self.sampling_matrix = np.array([
            [float(param_x[1]), float(param_x[2]), float(param_x[3])],
            [float(param_y[1]), float(param_y[2]), float(param_y[3])],
            [float(param_z[1]), float(param_z[2]), float(param_z[3])]
        ])
        
        self.sampling_step.append(float(param_x[1]))
        self.sampling_step.append(float(param_y[2]))
        self.sampling_step.append(float(param_z[3]))
        
    def read_atom_coordinates(self, cube_file):
        for i in range(6, self.num_atoms+6):
            element = dict_atom_number[int(util.remove_space_keep_num(self._lines[i])[0])]
            position = [float(pos) for pos in util.remove_space_keep_num(self._lines[i])[2:]]
            self.atom_coordinates[f'{element}{i-6}'] = position
            
    def find_system_origin(self):
        return 0

    def read_volume_data(self, cube_file):
        for i in range(self.num_atoms+6, len(self._lines)):
            data_i = util.remove_space_keep_num(self._lines[i])
            for data_ij in data_i:
                print(data_ij)

    def get_num_atoms(self):
        return self.num_atoms

    def get_sampling_step(self):
        return self.sampling_step

    def get_origin(self):
        return self.origin