import os
import sys

sys.path.insert(1, './bagelfit/src')
sys.path.insert(1, './examples/')


current_file_path = os.path.abspath(__file__)
project_root = os.path.dirname(os.path.dirname(os.path.dirname(current_file_path)))
sys.path.append(os.path.join(project_root,'./src/'))
sys.path.append(os.path.join(project_root,'./examples/'))


from BagelFitter import BagelFitter

from Torus import Torus
