
import IMP
import IMP.em
from sys import exit
from itertools import product, islice
from numpy import arange, sqrt
from time import time
import matplotlib.pyplot as plt
from Torus import Torus
import sys

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
	def __init__(self):
		"""
		Initializes the BagelFitter with necessary attributes.
		"""
		self.voxel_size = None
		self.dmap_out = None
		self.dmap_out_binary_flag = None

	def load_exprimental_map(self, input_map_path: str, voxel_size: int | None = None) -> None:
		"""
		Loads an experimental density map for processing.
		
		Args:
			input_map_path (str): Path to the input density map file.
			voxel_size (int, optional): Size of each voxel (default is determined from the map).
		"""
		self.input_map_path = input_map_path
		
		try:
			self.dmap = IMP.em.read_map(self.input_map_path)
			
			print("Input Map path:",self.input_map_path)

		except Exception as e:
			print(e)
			exit(1)

		if voxel_size==None:
			self.voxel_size = self.dmap.get_spacing()
		else:
			self.voxel_size = voxel_size
		print("Voxel size of the map:",self.voxel_size)

	def create_blank_density_map(self, n_voxels: int) -> IMP.em.DensityMap:
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

	def create_rectangular_blank_density_map(self, rec_bb=[(-100,-100,-100),(100,100,100)]) -> IMP.em.DensityMap:
		"""
		Creates an empty density map centered at (0,0,0).
		
		Args:
			n_voxels (int): Number of voxels in each dimension.
		
		Returns:
			IMP.em.DensityMap: A blank density map with specified voxel size.
		"""
		##l_low, l_high = -int(n_voxels/2.) * self.voxel_size, int(n_voxels/2.) * self.voxel_size

		lx_low, ly_low, lz_low = rec_bb[0][0], rec_bb[0][1], rec_bb[0][2]
		lx_high, ly_high, lz_high = rec_bb[1][0], rec_bb[1][1], rec_bb[1][2]
		bb_new = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(lx_low, ly_low, lz_low), IMP.algebra.Vector3D(lx_high, ly_high, lz_high))
		return IMP.em.create_density_map(bb_new, self.voxel_size)

	def calculate_dice_coefficient(self, dmap1: IMP.em.DensityMap, dmap2: IMP.em.DensityMap) -> float:
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

	def plot_voxel_values(self) -> None:
		"""
		Plots the histogram of voxel intensity values in the density map.
		"""
		vals = [self.dmap.get_value(vox) for vox in range(self.dmap.get_header().get_number_of_voxels())]
		plt.hist(vals, log=True)
		plt.xlabel("Value")
		plt.ylabel("Frequency (log scale)")
		plt.show()

	def fill_binary_density(self, torus: Torus) -> None:
		"""
		Generates a binary density map based on the torus parameters.
		
		Args:
			torus (Torus): Torus object used to define the density map.
		"""

		num_vox = self.dmap_out.get_header().get_number_of_voxels()

		xy_max_rad = torus.R + torus.extension
		print(torus.extension)

		for vox in range(0,num_vox):

			#Getting individual values per dimension
			x = self.dmap_out.get_location_in_dim_by_voxel(vox, 0) + self.voxel_size/2
			y = self.dmap_out.get_location_in_dim_by_voxel(vox, 1) + self.voxel_size/2
			z = self.dmap_out.get_location_in_dim_by_voxel(vox, 2) + self.voxel_size/2

			if torus.contains_point(x, y, z):
				cx, cy, _ = torus.center
				if sqrt((x - cx)**2 + (y - cy)**2) < xy_max_rad:
				#if sqrt(x**2 + y**2) < xy_max_rad:
					self.dmap_out.set_value(vox, 1.0)
				else:
					self.dmap_out.set_value(vox, 0.0)
			else:
				self.dmap_out.set_value(vox, 0.0)

	def fill_binary_density_with_multiple_torus_OLD(self, Tori: list) -> None:

		num_vox = self.dmap_out.get_header().get_number_of_voxels()

		for vox in range(0,num_vox):

			#Getting individual values per dimension
			x = self.dmap_out.get_location_in_dim_by_voxel(vox, 0) + self.voxel_size/2
			y = self.dmap_out.get_location_in_dim_by_voxel(vox, 1) + self.voxel_size/2
			z = self.dmap_out.get_location_in_dim_by_voxel(vox, 2) + self.voxel_size/2

			inside_any = False
			for torus in Tori:

				xy_max_rad = torus.R + torus.extension
				cx, cy, _ = torus.center

				if torus.contains_point(x, y, z):
					
					if sqrt((x - cx)**2 + (y - cy)**2) < xy_max_rad:
						inside_any = True
						break;

			self.dmap_out.set_value(vox, 1.0 if inside_any else 0.0)

	def fill_binary_density_with_multiple_torus(self, Tori: list) -> None:

		num_vox = self.dmap_out.get_header().get_number_of_voxels()

		for vox in range(num_vox):
			x = self.dmap_out.get_location_in_dim_by_voxel(vox, 0) + self.voxel_size/2
			y = self.dmap_out.get_location_in_dim_by_voxel(vox, 1) + self.voxel_size/2
			z = self.dmap_out.get_location_in_dim_by_voxel(vox, 2) + self.voxel_size/2

			# pick nearest torus in XY (by rho)
			best_torus = None
			best_rho = float("inf")

			for torus in Tori:
				cx, cy, _ = torus.center
				rho = sqrt((x - cx)**2 + (y - cy)**2)
				if rho < best_rho:
					best_rho = rho
					best_torus = torus

			# now evaluate only that torus
			xy_max_rad = best_torus.R + best_torus.extension

			inside = False
			if best_rho <= xy_max_rad and best_torus.contains_point(x, y, z):
				inside = True

			self.dmap_out.set_value(vox, 1.0 if inside else 0.0)

	

	def fill_nonbinary_density(self, torus: Torus) -> None:
		"""
		Generates a non-binary density (extracted from the input map) map based on the torus parameters.
		
		Args:
			torus (Torus): Torus object used to define the density map.
		"""

		try:
			num_vox = self.dmap.get_header().get_number_of_voxels()
		except Exception as e:
			print(e)
			print("Input map is required for filling torus with non-binary values")


		xy_max_rad = torus.R + torus.extension

		for vox in range(0,num_vox):

			#Getting individual values per dimension
			x = self.dmap.get_location_in_dim_by_voxel(vox, 0) + self.voxel_size/2
			y = self.dmap.get_location_in_dim_by_voxel(vox, 1) + self.voxel_size/2
			z = self.dmap.get_location_in_dim_by_voxel(vox, 2) + self.voxel_size/2

			if torus.contains_point(x, y, z):
				if sqrt(x**2 + y**2) < xy_max_rad:
					self.dmap_out.set_value(vox, self.dmap.get_value(vox))
				else:
					self.dmap_out.set_value(vox, 0.0)
			else:
				self.dmap_out.set_value(vox, 0.0)

	def fit_binary_torus(self, tor_R_range: tuple = (670, 680, 10), tor_r_range: tuple = (160, 180, 20), tor_th_range: tuple = (85.0, 95, 10), extension: float = 0.0) -> Torus:
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

		#This makes the datatype as python int rather than numpy int64. Increases the speed 50%
		param_sets = [(int(R), int(r), int(th)) for R in arange(*tor_R_range) for r in arange(*tor_r_range) for th in arange(*tor_th_range)]
		
		self.best_torus, max_cc_coef = None, -1
		
		start_time1 = time()

		for tor_R, tor_r, tor_th in param_sets:

			torus = Torus(tor_R, tor_r, tor_th, extension)
			self.fill_binary_density(torus)
			
			cc_coef = IMP.em.bayesem3d_get_cross_correlation_coefficient(self.dmap, self.dmap_out)
			
			if cc_coef > max_cc_coef:
				max_cc_coef = cc_coef
				self.best_torus = torus
				print("Best CC:", cc_coef, (tor_R, tor_r, tor_th))

		print(f"Execution Time: {time() - start_time1:.4f} seconds")

		if self.best_torus:
			print(f"Best Parameters: R={self.best_torus.R}, r={self.best_torus.r}, thickness={self.best_torus.thickness}")
			print(f"Highest CC Coefficient: {max_cc_coef:.4f}")
			self.best_torus.dmap = self.dmap_out
			return self.best_torus

	def fit_nonbinary_torus(self, tor_R_range: tuple = (670, 680, 10), tor_r_range: tuple = (160, 180, 20), tor_th_range: tuple = (85.0, 95, 10), extension: float = 0.0) -> Torus:
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

	def generate_binary_torus(self, tor_R: float, tor_r: float, tor_th: float, extension: float = 0.0, boundingbox_length: float = 2240, voxel_size: float = 10, outmap_fname: str = "torus_yeast_fitted.mrc", writemap: bool = True) -> Torus:
		"""
		Generates and writes a binary torus density map based on input parameters.
		
		Args:
			tor_R (float): Major radius of the torus.
			tor_r (float): Minor radius of the torus.
			tor_th (float): Thickness of the bilipid layer.
			extension (float, optional): Additional extension factor (default is 0.0).
			boundingbox_length (float, optional): Length of the bounding box for Torus centered at (0,0,0). Default value 2240 Å.
			voxel_size (float, optional): Individual voxel size in the output map file. Default value 10 Å.
			outmap_fname (str, optional): Output map file name of the torus. Default is "torus_yeast_fitted.mrc" in the current directory.
		
		Returns:
			Torus: Torus object based on input parameters.
		"""
		self.dmap_out_binary_flag = True

		self.voxel_size = voxel_size
		self.bb_length = boundingbox_length

		self.n_voxels = int(self.bb_length / self.voxel_size)

		self.dmap_out = self.create_blank_density_map(self.n_voxels)

		self.best_torus= None

		start_time2 = time()
		for tor_R, tor_r, tor_th in [(tor_R, tor_r, tor_th)]:

			torus = Torus(tor_R, tor_r, tor_th, extension)
			self.fill_binary_density(torus)

		self.best_torus = torus

		if writemap==True:
			IMP.em.write_map(self.dmap_out, outmap_fname)
		
		print(f"Torus with Parameters: R={self.best_torus.R}, r={self.best_torus.r}, thickness={self.best_torus.thickness}")
		print(f"Execution Time: {time() - start_time2:.4f} seconds")
		
		return self.best_torus

	def generate_multiple_binary_torus(self, tor_R: float, tor_r: float, tor_th: float, extension: float = 0.0, rectangular_bb = [(-150,-150,-150),(150,150,150)], torus_centers: list = [(0.0, 0.0, 0.0)], voxel_size: float = 10, outmap_fname: str = "torus_yeast_fitted.mrc", writemap: bool = True) -> Torus:
		"""
		Generates and writes a binary torus density map based on input parameters.
		
		Args:
			tor_R (float): Major radius of the torus.
			tor_r (float): Minor radius of the torus.
			tor_th (float): Thickness of the bilipid layer.
			extension (float, optional): Additional extension factor (default is 0.0).
			boundingbox_length (float, optional): Length of the bounding box for Torus centered at (0,0,0). Default value 2240 Å.
			voxel_size (float, optional): Individual voxel size in the output map file. Default value 10 Å.
			outmap_fname (str, optional): Output map file name of the torus. Default is "torus_yeast_fitted.mrc" in the current directory.
		
		Returns:
			Torus: Torus object based on input parameters.
		"""
		self.dmap_out_binary_flag = True

		self.voxel_size = voxel_size
		
		#rectangular_bb=[(-1500,-1500,-1500),(1500,1500,1500)]

		#self.bb_length = boundingbox_length
		#self.n_voxels = int(self.bb_length / self.voxel_size)

		#self.dmap_out = self.create_blank_density_map(self.n_voxels)
		self.dmap_out = self.create_rectangular_blank_density_map(rectangular_bb)

		self.best_torus= None

		start_time2 = time()

		TORI = []
		for tori_center in torus_centers:
			Tori_i = Torus(tor_R, tor_r, tor_th, extension, center=tori_center)
			TORI.append(Tori_i)

		self.fill_binary_density_with_multiple_torus(TORI)

		#self.best_torus = torus

		if writemap==True:
			IMP.em.write_map(self.dmap_out, outmap_fname)
		
		#print(f"Torus with Parameters: R={self.best_torus.R}, r={self.best_torus.r}, thickness={self.best_torus.thickness}")
		print(f"Execution Time: {time() - start_time2:.4f} seconds")
		
		return self.best_torus

	def generate_nonbinary_torus(self, tor_R: float, tor_r: float, tor_th: float, extension: float = 0.0, boundingbox_length: float = 2240, voxel_size: float = 10, outmap_fname: str = "torus_yeast_fitted.mrc", writemap: bool = True) -> Torus:
		"""
		Generates and writes a binary torus density map based on input parameters.
		
		Args:
			tor_R (float): Major radius of the torus.
			tor_r (float): Minor radius of the torus.
			tor_th (float): Thickness of the bilipid layer.
			extension (float, optional): Additional extension factor (default is 0.0).
			boundingbox_length (float, optional): Length of the bounding box for Torus centered at (0,0,0). Default value 2240 Å.
			voxel_size (float, optional): Individual voxel size in the output map file. Default value 10 Å.
			outmap_fname (str, optional): Output map file name of the torus. Default is "torus_yeast_fitted.mrc" in the current directory.
		
		Returns:
			Torus: Torus object based on input parameters.
		"""
		self.dmap_out_binary_flag = False

		self.voxel_size = voxel_size
		self.bb_length = boundingbox_length

		self.n_voxels = int(self.bb_length / self.voxel_size)

		self.dmap_out = self.create_blank_density_map(self.n_voxels)

		self.best_torus= None

		start_time2 = time()
		for tor_R, tor_r, tor_th in [(tor_R, tor_r, tor_th)]:

			torus = Torus(tor_R, tor_r, tor_th, extension)
			self.fill_nonbinary_density(torus)

		self.best_torus = torus

		if writemap==True:
			IMP.em.write_map(self.dmap_out, outmap_fname)
		
		print(f"Torus with Parameters: R={self.best_torus.R}, r={self.best_torus.r}, thickness={self.best_torus.thickness}")
		print(f"Execution Time: {time() - start_time2:.4f} seconds")
		
		return self.best_torus

	def write_torusmap_to_file(self, outmap_fname: str) -> None:
		"""
		Saves the torus density map to a file.
		
		Args:
			outmap_fname (str): Filename to save the torus density map.
		"""
		if self.best_torus==None:
			print("Best Torus object is None! ")
			exit(1)

		print(f"Saving torus map in file:{outmap_fname}")
		if self.dmap_out_binary_flag == True:
			self.fill_binary_density(self.best_torus)
		else:
			self.fill_nonbinary_density(self.best_torus)
		IMP.em.write_map(self.dmap_out, outmap_fname)
		return

	def score_torus_maps(self, map1: str, map2: str) -> float:
		"""
		Computes the cross-correlation coefficient between two torus maps.
		
		Args:
			map1 (str): Path to the first torus density map.
			map2 (str): Path to the second torus density map.
		
		Returns:
			float: Cross-correlation coefficient indicating similarity between the two maps.
		"""
		self.dmap1 = IMP.em.read_map(map1)

		self.dmap2 = IMP.em.read_map(map2)

		cc_coef = IMP.em.bayesem3d_get_cross_correlation_coefficient(self.dmap1, self.dmap2)

		print("Bayesian 3D crosscorrelation coeeficient",cc_coef)
		
		return cc_coef





