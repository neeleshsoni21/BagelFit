

from BagelFitter import BagelFitter

if __name__ == '__main__':
	
	#Below map has 11 million voxels    
	fitter= BagelFitter("./yeast_membrane/Yeast_C8_Double_MR_center.mrc")

	print("Map path:",fitter.input_map_path)
	print("Voxel size of the map:",fitter.voxel_size)

	tor_R_range=(670, 680, 10); tor_r_range=(160, 180, 20); tor_th_range=(85.0, 95, 10)
	extension = 0.0
	
	best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)

	#best_torus = fitter.fit_nonbinary_torus(tor_R_range, tor_r_range, tor_th_range, extension)

	fitter.write_torusmap_to_file("./yeast_membrane/torus_yeast_fitted.mrc")






