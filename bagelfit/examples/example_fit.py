
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
bagelfit_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(bagelfit_root)

#from bagelfit.src.BagelFitter import BagelFitter
#from bagelfit.src.Torus import Torus

import os
import bagelfit as bf

data_path = os.path.join(script_dir,"./yeast_membrane/")
#-----------------------------------------------------#
# Write a torus map using Torus parameters
#-----------------------------------------------------#
fitter= bf.BagelFitter()

tor_R = 670; tor_r=160; tor_th=85; extension=0.0
best_torus = fitter.generate_binary_torus(
	tor_R, tor_r, tor_th, 
	extension=0.0, 
	boundingbox_length=2240, 
	voxel_size=10.0, 
	outmap_fname=os.path.join(data_path,"torus_yeast_fitted.mrc" ))

#-----------------------------------------------------#
# Score a torus map with other torus maps or experimental maps
#-----------------------------------------------------#
fitter= bf.BagelFitter()

mapfile1 = os.path.join(data_path,"Yeast_C8_Double_MR_center.mrc" )
mapfile2 = os.path.join(data_path,"torus_yeast_fitted.mrc" )
fitter.score_torus_maps(mapfile1, mapfile2)


#-----------------------------------------------------#
# EXAMPLE Fit several torus onto nuclear membrane 
# and write the highes overlapping torus into the file
#-----------------------------------------------------#

fitter= bf.BagelFitter()

#Below map has 11 million voxels
fitter.load_exprimental_map(os.path.join(data_path,"Yeast_C8_Double_MR_center.mrc" ))

tor_R_range=(670, 680, 10); tor_r_range=(160, 180, 20); tor_th_range=(85, 95, 10)
extension = 0.0

best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)
#best_torus = fitter.fit_nonbinary_torus(tor_R_range, tor_r_range, tor_th_range, extension)

fitter.write_torusmap_to_file(os.path.join(data_path,"torus_yeast_fitted.mrc" ))











