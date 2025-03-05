
import IMP
import IMP.em
from itertools import product
from numpy import arange, sqrt
from time import time
import matplotlib.pyplot as plt
from Torus import Torus

class BagelFitter:
	"""
	Fits a torus onto a nuclear membrane by searching for the best parameters.
	
	Attributes:
		best_torus (Torus): The best-fitting torus found during the fitting process.
		dmap (IMP.em.DensityMap): Input density map of the nuclear membrane.
		dmap_out (IMP.em.DensityMap): Output density map after fitting.
		dmap_out_binary_flag (bool): Flag indicating if the density map is binary.
		input_map_path (str): Path to the input density map file.
		voxel_size (int): Size of each voxel in the density map (default is 10).
	"""
	def __init__(self, input_map_path, voxel_size=10):
		"""
		Initializes the NuTorusFitter with an input density map.
		
		Args:
			input_map_path (str): Path to the input density map file.
			voxel_size (int, optional): Size of each voxel in the density map (default is 10).
		"""
		self.input_map_path = input_map_path
		self.voxel_size = voxel_size
		self.dmap = IMP.em.read_map(input_map_path)
		self.dmap_out = None
		self.dmap_out_binary_flag = None


	def create_blank_density_map(self, n_voxels):
		"""
		Creates an empty density map centered at (0,0,0).
		
		Args:
			n_voxels (int): Number of voxels in each dimension.
		
		Returns:
			IMP.em.DensityMap: A blank density map with specified voxel size.
		"""
		l_low, l_high = -int(n_voxels/2.) * self.voxel_size, int(n_voxels/2.) * self.voxel_size
		bb_new = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(l_low, l_low, l_low), IMP.algebra.Vector3D(l_high, l_high, l_high))
		return IMP.em.create_density_map(bb_new, self.voxel_size)

	def calculate_dice_coefficient(self, dmap1, dmap2):
		"""
		Computes the Dice Coefficient (F1 Score) for binary overlap between two density maps.
		
		Args:
			dmap1 (IMP.em.DensityMap): First density map.
			dmap2 (IMP.em.DensityMap): Second density map.
		
		Returns:
			float: Dice Coefficient representing the overlap of the two maps.
		"""
		num_voxels = dmap1.get_header().get_number_of_voxels()
		I_dmap1, I_dmap2, I_dmap1_dmap2 = 0, 0, 0

		for vox in range(num_voxels):
			dmap1_val, dmap2_val = dmap1.get_value(vox), dmap2.get_value(vox)
			I_dmap1 += dmap1_val == 1
			I_dmap2 += dmap2_val == 1
			I_dmap1_dmap2 += dmap1_val * dmap2_val == 1

		return 2 * I_dmap1_dmap2 / (I_dmap1 + I_dmap2)

	def plot_voxel_values(self):
		"""
		Plots the histogram of voxel intensity values in the density map.
		"""
		vals = [self.dmap.get_value(vox) for vox in range(self.dmap.get_header().get_number_of_voxels())]
		plt.hist(vals, log=True)
		plt.xlabel("Value")
		plt.ylabel("Frequency (log scale)")
		plt.show()

	def fill_binary_density(self, torus):
		"""
		Generates a binary density map based on the torus parameters.
		
		Args:
			torus (Torus): Torus object used to define the density map.
		"""

		num_vox = self.dmap.get_header().get_number_of_voxels()
		#num_vox = 3000000

		xy_max_rad = torus.R + torus.extension

		for vox in range(0,num_vox):

			#Getting individual values per dimension
			x = self.dmap.get_location_in_dim_by_voxel(vox, 0)
			y = self.dmap.get_location_in_dim_by_voxel(vox, 1)
			z = self.dmap.get_location_in_dim_by_voxel(vox, 2)

			if torus.contains_point(x, y, z):
				if sqrt(x**2 + y**2) < xy_max_rad:
					self.dmap_out.set_value(vox, 1.0)
				else:
					self.dmap_out.set_value(vox, 0.0)
			else:
				self.dmap_out.set_value(vox, 0.0)

	def fill_nonbinary_density(self, torus):
		"""
		Generates a non-binary density map based on the torus parameters.
		
		Args:
			torus (Torus): Torus object used to define the density map.
		"""

		num_vox = self.dmap.get_header().get_number_of_voxels()

		xy_max_rad = torus.R + torus.extension

		for vox in range(0,num_vox):

			#Getting individual values per dimension
			x = self.dmap.get_location_in_dim_by_voxel(vox, 0)
			y = self.dmap.get_location_in_dim_by_voxel(vox, 1)
			z = self.dmap.get_location_in_dim_by_voxel(vox, 2)

			if torus.contains_point(x, y, z):
				if sqrt(x**2 + y**2) < xy_max_rad:
					self.dmap_out.set_value(vox, self.dmap.get_value(vox))
				else:
					self.dmap_out.set_value(vox, 0.0)
			else:
				self.dmap_out.set_value(vox, 0.0)

	def fit_binary_torus(self, tor_R_range=(670, 680, 10), tor_r_range=(160, 180, 20), tor_th_range=(85.0, 95, 10), extension=0.0):
		"""
		Fits a binary torus to the density map by searching for optimal parameters.
		
		Args:
			tor_R_range (tuple, optional): Range of major radius values (start, stop, step).
			tor_r_range (tuple, optional): Range of minor radius values (start, stop, step).
			tor_th_range (tuple, optional): Range of thickness values (start, stop, step).
			extension (float, optional): Additional extension factor (default is 0.0).
		
		Returns:
			Torus: Best-fitting torus object based on maximum cross-correlation coefficient.
		"""
		self.dmap_out_binary_flag = True

		bb = IMP.em.get_bounding_box(self.dmap)
		n_voxels = (bb[1][0] - bb[0][0]) / self.voxel_size
		self.dmap_out = self.create_blank_density_map(n_voxels)

		param_sets = list(product(arange(*tor_R_range), arange(*tor_r_range), arange(*tor_th_range)))
		self.best_torus, max_cc_coef = None, -1

		for tor_R, tor_r, tor_th in param_sets:

			start_time = time()
			torus = Torus(tor_R, tor_r, tor_th, extension)
			self.fill_binary_density(torus)

			cc_coef = IMP.em.bayesem3d_get_cross_correlation_coefficient(self.dmap, self.dmap_out)
			if cc_coef > max_cc_coef:
				max_cc_coef = cc_coef
				self.best_torus = torus
				print("Best CC:", cc_coef, (tor_R, tor_r, tor_th))

			print(f"Execution Time: {time() - start_time:.4f} seconds")

		if self.best_torus:
			print(f"Best Parameters: R={self.best_torus.R}, r={self.best_torus.r}, thickness={self.best_torus.thickness}")
			print(f"Highest CC Coefficient: {max_cc_coef:.4f}")
			self.best_torus.dmap = self.dmap_out
			
			return self.best_torus

	def fit_nonbinary_torus(self, tor_R_range=(670, 680, 10), tor_r_range=(160, 180, 20), tor_th_range=(85.0, 95, 10), extension=0.0):
		"""
		Fits a non-binary torus to the density map by searching for optimal parameters.
		
		Args:
			tor_R_range (tuple, optional): Range of major radius values (start, stop, step).
			tor_r_range (tuple, optional): Range of minor radius values (start, stop, step).
			tor_th_range (tuple, optional): Range of thickness values (start, stop, step).
			extension (float, optional): Additional extension factor (default is 0.0).
		
		Returns:
			Torus: Best-fitting torus object based on maximum cross-correlation coefficient.
		"""
		self.dmap_out_binary_flag = False

		bb = IMP.em.get_bounding_box(self.dmap)
		n_voxels = (bb[1][0] - bb[0][0]) / self.voxel_size
		self.dmap_out = self.create_blank_density_map(n_voxels)

		param_sets = list(product(arange(*tor_R_range), arange(*tor_r_range), arange(*tor_th_range)))
		self.best_torus, max_cc_coef = None, -1

		for tor_R, tor_r, tor_th in param_sets:

			start_time = time()
			torus = Torus(tor_R, tor_r, tor_th, extension)
			self.fill_nonbinary_density(torus)

			cc_coef = IMP.em.bayesem3d_get_cross_correlation_coefficient(self.dmap, self.dmap_out)
			if cc_coef > max_cc_coef:
				max_cc_coef = cc_coef
				self.best_torus = torus
				print("Best CC:", cc_coef, (tor_R, tor_r, tor_th))

			print(f"Execution Time: {time() - start_time:.4f} seconds")

		if self.best_torus:
			print(f"Best Parameters: R={self.best_torus.R}, r={self.best_torus.r}, thickness={self.best_torus.thickness}")
			print(f"Highest CC Coefficient: {max_cc_coef:.4f}")
			self.best_torus.dmap = self.dmap_out
			
			return self.best_torus

	def write_torusmap_to_file(self, outmap_fname):
		"""
		Saves the torus density map to a file.
		
		Args:
			outmap_fname (str): Filename to save the torus density map.
		"""
		print(f"Saving torus map in file:{outmap_fname}")
		if self.dmap_out_binary_flag == True:
			self.fill_binary_density(self.best_torus)
		else:
			self.fill_nonbinary_density(self.best_torus)
		IMP.em.write_map(self.dmap_out, outmap_fname)
		return
