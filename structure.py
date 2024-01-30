import numpy as np

from grades.const import dict_atom_number
import util

class Structure():
    def __init__(self, file_path):
        self.file_path = file_path
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
        
    def parse_cube_file(self):
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