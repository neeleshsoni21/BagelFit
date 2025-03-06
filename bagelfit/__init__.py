import os
import sys

#sys.path.insert(1, './bagelfit/src')
#sys.path.insert(1, './examples/')


current_file_path = os.path.abspath(__file__)
#project_root = os.path.dirname(os.path.dirname(os.path.dirname(current_file_path)))
project_root = os.path.dirname(os.path.dirname(current_file_path))
sys.path.append(os.path.join(project_root,'./bagelfit/src/'))
sys.path.append(os.path.join(project_root,'./bagelfit/examples/'))


from BagelFitter import BagelFitter

from Torus import Torus
