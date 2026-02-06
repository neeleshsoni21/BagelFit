"""
Example script for generating and scoring torus maps in nuclear membrane fitting.

"""
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
bagelfit_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(bagelfit_root)

def generate_binary_torus_map(boundingbox, extension):
	"""
	Generates a binary torus map using predefined torus parameters and saves it to a file.
	
	Returns:
		None

	Examples:
		.. code-block:: python

			data_path = os.path.join(script_dir,"./yeast_membrane/")

			fitter= bf.BagelFitter()

			tor_R = 660; tor_r=140; tor_th=55; extension=0.0

			best_torus = fitter.generate_binary_torus(
				tor_R, tor_r, tor_th, 
				extension=0.0, 
				boundingbox_length=2240, 
				voxel_size=10.0, 
				outmap_fname=os.path.join(data_path,"torus_yeast_fitted.mrc" ))
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir,"./yeast_membrane/")
	print(data_path)
	#-----------------------------------------------------#
	# Write a torus map using Torus parameters
	#-----------------------------------------------------#
	fitter= bf.BagelFitter()

	tor_R = 650; tor_r=125; tor_th=40;
	best_torus = fitter.generate_binary_torus(
		tor_R, tor_r, tor_th, 
		extension=extension, 
		boundingbox_length=boundingbox, 
		voxel_size=10.0, 
		outmap_fname=os.path.join(data_path,"torus_yeast_fitted.mrc" ))

	return

def generate_two_binary_torus_map():

	import bagelfit as bf

	data_path = os.path.join(script_dir,"./yeast_membrane/")
	print(data_path)
	#-----------------------------------------------------#
	# Write a torus map using Torus parameters
	#-----------------------------------------------------#
	fitter= bf.BagelFitter()

	tor_R = 650; tor_r=125; tor_th=40;
	Inter_NPC_distance = 1500
	ext_distance = Inter_NPC_distance/2.0 - tor_R
	boundingbox = 3000
	#THIS extension to torus will cover the entire rectangular bounding box of length: Inter_NPC_distance
	#ext_R = (2*(tor_R+ext_distance)**2)**0.5 
	ext_R = 1000.0
	
	best_torus = fitter.generate_multiple_binary_torus(
		tor_R, tor_r, tor_th, 
		extension=ext_R, 
		rectangular_bb = [(-1500,-750,-1000),(1500,750,500)],
		#torus_centers=[(-1500, 0.0, 0.0),(0.0, 0.0, 0.0),(1500, 0.0, 0.0)],
		torus_centers=[(-750.0, 0.0, 0.0),(750, 0.0, 0.0)],
		voxel_size=10.0, 
		outmap_fname=os.path.join(data_path,"multi_torus_yeast.mrc" ))

	return

def generate_four_binary_torus_map():

	import bagelfit as bf

	data_path = os.path.join(script_dir,"./yeast_membrane/")
	print(data_path)
	#-----------------------------------------------------#
	# Write a torus map using Torus parameters
	#-----------------------------------------------------#
	fitter= bf.BagelFitter()

	tor_R = 650; tor_r=125; tor_th=40;
	#Inter_NPC_distance = 1500
	#ext_distance = Inter_NPC_distance/2.0 - tor_R

	#THIS extension to torus will cover the entire rectangular bounding box of length: Inter_NPC_distance
	#ext_R = (2*(tor_R+ext_distance)**2)**0.5 
	ext_R = 1000.0
	
	best_torus = fitter.generate_multiple_binary_torus(
		tor_R, tor_r, tor_th, 
		extension=ext_R, 
		rectangular_bb = [(-1500,-2200,-1000),(1500,2200,500)],
		#torus_centers=[(-1500, 0.0, 0.0),(0.0, 0.0, 0.0),(1500, 0.0, 0.0)],
		torus_centers=[(-750.0, 0.0, 0.0),(750, 0.0, 0.0),(0.0, -1300.0, 0.0),(0.0, 1300, 0.0)],
		voxel_size=10.0, 
		outmap_fname=os.path.join(data_path,"four_torus_yeast.mrc" ))

	return

def generate_nonbinary_torus_map(experimental_map, boundingbox, extension):
	"""
	Generates a binary torus map using predefined torus parameters and saves it to a file.
	
	Returns:
		None

	Examples:
		.. code-block:: python

			data_path = os.path.join(script_dir,"./yeast_membrane/")

			fitter= bf.BagelFitter()

			tor_R = 660; tor_r=140; tor_th=55; extension=0.0

			best_torus = fitter.generate_nonbinary_torus(
				tor_R, tor_r, tor_th, 
				extension=0.0, 
				boundingbox_length=2240, 
				voxel_size=10.0, 
				outmap_fname=os.path.join(data_path,"torus_yeast_fitted.mrc" ))
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir,"./yeast_membrane/")
	print(data_path)
	#-----------------------------------------------------#
	# Write a torus map using Torus parameters
	#-----------------------------------------------------#
	fitter= bf.BagelFitter()

	fitter.load_exprimental_map(
		input_map_path=experimental_map,
		voxel_size=10.0)

	tor_R = 650; tor_r=125; tor_th=40;
	best_torus = fitter.generate_nonbinary_torus(
		tor_R, tor_r, tor_th, 
		extension=extension, 
		boundingbox_length=boundingbox, 
		voxel_size=10.0, 
		outmap_fname=os.path.join(data_path,"torus_yeast_fitted.mrc" ))

	return


	generate_nonbinary_torus

def score_torus_map_with_experimental_map():
	"""
	Compares a generated torus map with an experimental map using a scoring function.
	
	Returns:
		None

	Examples:
		.. code-block:: python

			fitter= bf.BagelFitter()

			mapfile1 = os.path.join(data_path,"Yeast_C8_Double_MR_center.mrc" )
			
			mapfile2 = os.path.join(data_path,"torus_yeast_fitted.mrc" )
			
			fitter.score_torus_maps(mapfile1, mapfile2)

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

	Examples:
		.. code-block:: python

			fitter= bf.BagelFitter()
			
			fitter.load_exprimental_map(os.path.join(data_path,"Yeast_C8_Double_MR_center.mrc" ))
			
			tor_R_range=(660, 670, 10); 
			tor_r_range=(140, 160, 20); 
			tor_th_range=(55, 65, 10); 
			extension = 0.0
			
			best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)
			
			fitter.write_torusmap_to_file(os.path.join(data_path,"torus_yeast_fitted.mrc" ))
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

	#---------------------------------
	# EXAMPLE CODE to generate torus and score using experimental 
	# maps and bet the best toroid map and parameters
	#generate_binary_torus_map(boundingbox=2500, extension=1000)
	
	#generate_nonbinary_torus_map(
		#experimental_map='./yeast_membrane/Yeast_C8_Double_MR_center_resampled.mrc',
		#boundingbox=3000, extension=0) #Dimension cannot be more than the experimental map for non-binary points

	#score_torus_map_with_experimental_map()

	#generate_bestfit_torus_map()

	#---------------------------------


	#---------------------------------
	# EXAMPLE code to generate extended membrane around the torus
	#---------------------------------
	#generate_binary_torus_map(boundingbox=2500, extension=1768) #(2(1250**2))**0.5 = 1768


	#---------------------------------
	# EXAMPLE code to generate TWO torus with extended membrane around the torus
	#---------------------------------
	#generate_two_binary_torus_map()

	generate_four_binary_torus_map()


	











