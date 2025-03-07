"""
Exanple script for generating and scoring torus maps in nuclear membrane fitting.

"""
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
bagelfit_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(bagelfit_root)

def generate_binary_torus_map():
	"""
	Generates a binary torus map using predefined torus parameters and saves it to a file.
	
	Returns:
		None
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir,"./yeast_membrane/")
	#-----------------------------------------------------#
	# Write a torus map using Torus parameters
	#-----------------------------------------------------#
	fitter= bf.BagelFitter()


	tor_R = 660; tor_r=140; tor_th=55; extension=0.0
	best_torus = fitter.generate_binary_torus(
		tor_R, tor_r, tor_th, 
		extension=0.0, 
		boundingbox_length=2240, 
		voxel_size=10.0, 
		outmap_fname=os.path.join(data_path,"torus_yeast_fitted.mrc" ))

	return

def score_torus_map_with_experimental_map():
	"""
	Compares a generated torus map with an experimental map using a scoring function.
	
	Returns:
		None
	"""
	#-----------------------------------------------------#
	# Score a torus map with other torus maps or experimental maps
	#-----------------------------------------------------#
	fitter= bf.BagelFitter()

	mapfile1 = os.path.join(data_path,"Yeast_C8_Double_MR_center.mrc" )
	mapfile2 = os.path.join(data_path,"torus_yeast_fitted.mrc" )
	fitter.score_torus_maps(mapfile1, mapfile2)

	return


def generate_bestfit_torus_map():
	"""
	Fits several torus models onto the nuclear membrane and saves the best fit.
	
	Returns:
		None
	"""
	#-----------------------------------------------------#
	# EXAMPLE Fit several torus onto nuclear membrane 
	# and write the highes overlapping torus into the file
	#-----------------------------------------------------#

	fitter= bf.BagelFitter()

	#Below map has 11 million voxels
	fitter.load_exprimental_map(os.path.join(data_path,"Yeast_C8_Double_MR_center.mrc" ))

	tor_R_range=(660, 670, 10); tor_r_range=(140, 160, 20); tor_th_range=(55, 65, 10)
	extension = 0.0

	best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)
	#best_torus = fitter.fit_nonbinary_torus(tor_R_range, tor_r_range, tor_th_range, extension)

	fitter.write_torusmap_to_file(os.path.join(data_path,"torus_yeast_fitted.mrc" ))

	return

if __name__ == '__main__':

	generate_binary_torus_map()

	score_torus_map_with_experimental_map()

	generate_bestfit_torus_map()









