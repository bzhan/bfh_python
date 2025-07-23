""" Data for testing stored here """

import os

data_root_path = os.path.dirname(__file__)

def data_file(file_name):
    path = os.path.join(data_root_path, file_name)
    return open(path, "r")

def data_file_contents(file_name):
    file = data_file(file_name)
    ans = file.read()
    file.close()
    return ans

    
    
