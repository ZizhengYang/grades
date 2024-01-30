import re

def remove_space_keep_num(self, in_str):
    str_params = re.split(r"[;,\s]\s*", in_str)
    ret_array = []
    for str_param in str_params:
        if not str_param == '':
            ret_array.append(str_param)
    return ret_array