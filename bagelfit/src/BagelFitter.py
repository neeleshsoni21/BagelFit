
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
	Generate, score, and fit torus-based membrane models against an experimental density map.
	
	Core concepts:
	- A `Torus` defines a membrane *shell* (band of thickness `tor_th` within the tube radius `tor_r`).
	- *Binary* maps are occupancy masks (0/1) indicating whether each voxel is inside the shell.
	- *Non-binary* maps copy voxel densities from the experimental map inside the shell (0 outside).
	
	Fitting:
	- `fit_*_torus()` performs a grid search over (R, r, thickness) and selects the parameters
	  that maximize a similarity score (cross-correlation coefficient).
	"""
	def __init__(self):
		"""
		Initializes the BagelFitter with necessary attributes.
		"""
		self.voxel_size = None
		self.dmap_out = None
		self.dmap_out_binary_flag = None

	def load_exprimental_map(self, input_map_path, voxel_size= None):
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

	def create_rectangular_blank_density_map(self, rec_bb=[(-100,-100,-100),(100,100,100)]):
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
		Fill `self.dmap_out` with a *binary* torus membrane mask (0/1).
		
		A voxel is set to 1 when:
		- it lies inside the torus membrane shell (see `Torus.contains_point()`), and
		- its in-plane (XY) radius satisfies ρ_xy <= (R + extension). This cutoff is what
		  allows "extended membrane" coverage beyond the torus ring.
		
		Args:
		    torus: Torus defining the membrane shell and in-plane extension cutoff.
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

	def fill_binary_density_with_multiple_torus_OLD(self, Tori):

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

	def fill_binary_density_with_multiple_torus(self, Tori):
		"""
		Fill `self.dmap_out` with a binary mask for multiple tori.
		
		For each voxel we choose the nearest torus in the XY plane (by ρ_xy), then evaluate
		membership against that torus' membrane shell and its in-plane extension cutoff.
		This avoids double-counting when multiple tori overlap in space.
		
		Args:
		    Tori: List of `Torus` objects with potentially different centers.
		"""
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

	

	def fill_nonbinary_density(self, torus):
		"""
		Fill `self.dmap_out` with a *non-binary* torus map derived from the experimental map.
		
		For voxels inside the torus membrane shell (and within the in-plane extension cutoff),
		the output voxel value is copied from the experimental map. Voxels outside are set to 0.
		
		Args:
		    torus: Torus defining the membrane shell and in-plane extension cutoff.
		
		Raises:
		    RuntimeError: If no experimental map has been loaded.
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

	def fit_binary_torus(self, tor_R_range= (670, 680, 10), tor_r_range= (160, 180, 20), tor_th_range= (85.0, 95, 10), extension= 0.0):
		"""
		Fit a *binary* torus membrane model by grid-searching torus parameters.
		
		This generates candidate binary occupancy maps (0/1) and selects the parameters that
		maximize the cross-correlation coefficient against the experimental map.
		
		Args:
		    tor_R_range: (start, stop, step) for major radius R.
		    tor_r_range: (start, stop, step) for minor radius r (tube radius).
		    tor_th_range: (start, stop, step) for membrane-shell thickness.
		    extension: In-plane (XY) radial allowance used as ρ_xy <= (R + extension).
		
		Returns:
		    Best-fitting `Torus`.
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

	def fit_nonbinary_torus(self, tor_R_range= (670, 680, 10), tor_r_range= (160, 180, 20), tor_th_range= (85.0, 95, 10), extension= 0.0):
		"""
		Fit a *non-binary* torus membrane model by grid-searching torus parameters.
		
		This generates candidate non-binary maps (experimental densities inside the shell) and
		selects the parameters that maximize the cross-correlation coefficient.
		
		Args:
		    tor_R_range: (start, stop, step) for major radius R.
		    tor_r_range: (start, stop, step) for minor radius r (tube radius).
		    tor_th_range: (start, stop, step) for membrane-shell thickness.
		    extension: In-plane (XY) radial allowance used as ρ_xy <= (R + extension).
		
		Returns:
		    Best-fitting `Torus`.
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

	def generate_binary_torus(self, tor_R, tor_r, tor_th, extension= 0.0, boundingbox_length= 2240, voxel_size= 10, outmap_fname= "torus_yeast_fitted.mrc", writemap= True):
		"""
		Generate a *binary* torus occupancy map and optionally write it to disk.
		
		Args:
		    tor_R: Major radius R.
		    tor_r: Minor radius r (tube radius).
		    tor_th: Membrane-shell thickness within the tube.
		    extension: In-plane (XY) radial allowance used as ρ_xy <= (R + extension).
		    boundingbox_length: Side length of the (cubic) output bounding box (centered at origin).
		    voxel_size: Output voxel spacing.
		    outmap_fname: Output filename (.mrc).
		    writemap: If True, writes the map to `outmap_fname`.
		
		Returns:
		    Torus: The torus parameters used to generate the map.
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

	def generate_multiple_binary_torus(self, tor_R, tor_r, tor_th, extension= 0.0, rectangular_bb = [(-150,-150,-150),(150,150,150)], torus_centers= [(0.0, 0.0, 0.0)], voxel_size= 10, outmap_fname= "torus_yeast_fitted.mrc", writemap= True):
		"""
		Generate a *binary* occupancy map containing multiple tori and optionally write it to disk.
		
		Args:
		    tor_R: Major radius R (shared by all tori).
		    tor_r: Minor radius r (shared by all tori).
		    tor_th: Membrane-shell thickness (shared by all tori).
		    extension: In-plane (XY) radial allowance used as ρ_xy <= (R + extension) per torus.
		    rectangular_bb: Axis-aligned rectangular bounding box as [(xmin,ymin,zmin), (xmax,ymax,zmax)].
		    torus_centers: List of torus centers (cx, cy, cz) in map coordinates.
		    voxel_size: Output voxel spacing.
		    outmap_fname: Output filename (.mrc).
		    writemap: If True, writes the map to `outmap_fname`.
		
		Returns:
		    Torus: The *last* torus created (all share the same parameters but differ in center).
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

	def generate_nonbinary_torus(self, tor_R, tor_r, tor_th, extension= 0.0, boundingbox_length= 2240, voxel_size= 10, outmap_fname= "torus_yeast_fitted.mrc", writemap= True):
		"""
		Generate a *non-binary* torus map (experimental densities inside the shell) and optionally write it to disk.
		
		Requires that an experimental map has been loaded via `load_exprimental_map()`.
		
		Args:
		    tor_R: Major radius R.
		    tor_r: Minor radius r (tube radius).
		    tor_th: Membrane-shell thickness within the tube.
		    extension: In-plane (XY) radial allowance used as ρ_xy <= (R + extension).
		    boundingbox_length: Side length of the (cubic) output bounding box.
		    voxel_size: Output voxel spacing.
		    outmap_fname: Output filename (.mrc).
		    writemap: If True, writes the map to `outmap_fname`.
		
		Returns:
		    Torus: The torus parameters used to generate the map.
		"""		self.dmap_out_binary_flag = False

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

	def write_torusmap_to_file(self, outmap_fname):
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

	def score_torus_maps(self, map1, map2):
		"""
		Compute the cross-correlation coefficient between two density maps.
		
		Args:
		    map1: Path to the first map file (.mrc).
		    map2: Path to the second map file (.mrc).
		
		Returns:
		    Cross-correlation coefficient (float). Higher indicates greater similarity.
		"""
		self.dmap1 = IMP.em.read_map(map1)

		self.dmap2 = IMP.em.read_map(map2)

		cc_coef = IMP.em.bayesem3d_get_cross_correlation_coefficient(self.dmap1, self.dmap2)

		print("Bayesian 3D crosscorrelation coeeficient",cc_coef)
		
		return cc_coef





