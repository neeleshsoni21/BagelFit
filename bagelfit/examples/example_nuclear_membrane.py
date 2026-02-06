"""
Example script for generating and scoring torus maps in nuclear membrane fitting.

This script demonstrates how to:
- Generate single or multiple (binary) torus occupancy maps.
- Generate a non-binary (density/weighted) torus map aligned to an experimental map.
- Score a generated torus map against an experimental reference map.
- Fit torus parameters by grid search and export the best-fit torus map.

Notes:
- "Binary torus map" means voxels are labeled as inside/outside the torus (0/1).
- "Non-binary torus map" means voxels contain continuous density values (not just 0/1).
- `extension` is used to expand the modeled membrane region beyond the base torus.
"""
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
bagelfit_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(bagelfit_root)


def generate_binary_torus_map(boundingbox, extension):
	"""
	Generate a single *binary* torus occupancy map and save it to disk.

	Args:
		boundingbox (float | int):
			Side length of the (typically cubic) bounding box in the same units as
			the torus parameters (e.g., Angstroms).
		extension (float):
			Additional radial/planar extension used to expand the modeled membrane
			region around the torus (implementation-specific to bagelfit).

	Returns:
		None

	Example:
		.. code-block:: python

			data_path = os.path.join(script_dir, "./yeast_membrane/")
			fitter = bf.BagelFitter()

			tor_R = 660
			tor_r = 140
			tor_th = 55
			extension = 0.0

			best_torus = fitter.generate_binary_torus(
				tor_R, tor_r, tor_th,
				extension=extension,
				boundingbox_length=2240,
				voxel_size=10.0,
				outmap_fname=os.path.join(data_path, "torus_yeast_fitted.mrc"),
			)
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir, "./yeast_membrane/")
	print(data_path)

	# -----------------------------------------------------#
	# Generate and write a *binary* torus map from parameters
	# -----------------------------------------------------#
	fitter = bf.BagelFitter()

	tor_R = 650
	tor_r = 125
	tor_th = 40

	# Writes a binary occupancy map (.mrc) where voxels inside the torus region are 1.
	best_torus = fitter.generate_binary_torus(
		tor_R, tor_r, tor_th,
		extension=extension,
		boundingbox_length=boundingbox,
		voxel_size=10.0,
		outmap_fname=os.path.join(data_path, "torus_yeast_fitted.mrc"),
	)

	return


def generate_two_binary_torus_map():
	"""
	Generate a *binary* map containing two toroids (two centers) and save it to disk.

	This example uses a rectangular bounding box and places two torus centers along X.

	Returns:
		None
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir, "./yeast_membrane/")
	print(data_path)

	# -----------------------------------------------------#
	# Generate and write a *binary* multi-torus map
	# -----------------------------------------------------#
	fitter = bf.BagelFitter()

	tor_R = 650
	tor_r = 125
	tor_th = 40

	# Desired center-to-center spacing between toroids
	Inter_NPC_distance = 1500

	# How much extra "reach" is needed so the extended membrane spans the gap
	# between two toroids (conceptually: half-distance minus major radius).
	ext_distance = Inter_NPC_distance / 2.0 - tor_R

	# Overall map size control (not directly used below since rectangular_bb is provided).
	boundingbox = 3000

	# Extension used to expand the modeled membrane region around each torus.
	# This can be tuned to ensure coverage of the rectangular bounding box.
	# ext_R = (2*(tor_R+ext_distance)**2)**0.5  # optional derived estimate
	ext_R = 1000.0

	best_torus = fitter.generate_multiple_binary_torus(
		tor_R, tor_r, tor_th,
		extension=ext_R,
		# Explicit rectangular bounding box: [(xmin,ymin,zmin), (xmax,ymax,zmax)]
		rectangular_bb=[(-1500, -750, -1000), (1500, 750, 500)],
		# Two torus centers separated along X
		torus_centers=[(-750.0, 0.0, 0.0), (750.0, 0.0, 0.0)],
		voxel_size=10.0,
		outmap_fname=os.path.join(data_path, "multi_torus_yeast.mrc"),
	)

	return


def generate_four_binary_torus_map():
	"""
	Generate a *binary* map containing four toroids (four centers) and save it to disk.

	This example uses a rectangular bounding box and places torus centers in a cross-like
	configuration: two along X and two along Y.

	Returns:
		None
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir, "./yeast_membrane/")
	print(data_path)

	# -----------------------------------------------------#
	# Generate and write a *binary* multi-torus map
	# -----------------------------------------------------#
	fitter = bf.BagelFitter()

	tor_R = 650
	tor_r = 125
	tor_th = 40

	# Extension used to expand the modeled membrane region around each torus.
	# ext_R = (2*(tor_R+ext_distance)**2)**0.5  # optional derived estimate
	ext_R = 1000.0

	best_torus = fitter.generate_multiple_binary_torus(
		tor_R, tor_r, tor_th,
		extension=ext_R,
		# Explicit rectangular bounding box: [(xmin,ymin,zmin), (xmax,ymax,zmax)]
		rectangular_bb=[(-1500, -2200, -1000), (1500, 2200, 500)],
		# Four torus centers: two on X-axis and two on Y-axis
		torus_centers=[
			(-750.0, 0.0, 0.0),
			(750.0, 0.0, 0.0),
			(0.0, -1300.0, 0.0),
			(0.0, 1300.0, 0.0),
		],
		voxel_size=10.0,
		outmap_fname=os.path.join(data_path, "four_torus_yeast.mrc"),
	)

	return


def generate_nonbinary_torus_map(experimental_map, boundingbox, extension):
	"""
	Generate a *non-binary* (continuous-density) torus map and save it to disk.

	This function loads an experimental map first, then generates a non-binary torus
	map (e.g., weighted density rather than 0/1 occupancy).

	Args:
		experimental_map (str):
			Path to the experimental .mrc map used for alignment/resampling context.
		boundingbox (float | int):
			Side length of the output bounding box (must be compatible with the loaded
			experimental map and/or the library's constraints).
		extension (float):
			Additional extension used to expand the modeled membrane region around the torus.

	Returns:
		None

	Example:
		.. code-block:: python

			data_path = os.path.join(script_dir, "./yeast_membrane/")
			fitter = bf.BagelFitter()

			fitter.load_exprimental_map(
				input_map_path=experimental_map,
				voxel_size=10.0,
			)

			best_torus = fitter.generate_nonbinary_torus(
				tor_R, tor_r, tor_th,
				extension=extension,
				boundingbox_length=boundingbox,
				voxel_size=10.0,
				outmap_fname=os.path.join(data_path, "torus_yeast_fitted.mrc"),
			)
	"""
	import bagelfit as bf

	data_path = os.path.join(script_dir, "./yeast_membrane/")
	print(data_path)

	# -----------------------------------------------------#
	# Load experimental map, then generate a non-binary torus map
	# -----------------------------------------------------#
	fitter = bf.BagelFitter()

	fitter.load_exprimental_map(
		input_map_path=experimental_map,
		voxel_size=10.0,
	)

	tor_R = 650
	tor_r = 125
	tor_th = 40

	best_torus = fitter.generate_nonbinary_torus(
		tor_R, tor_r, tor_th,
		extension=extension,
		boundingbox_length=boundingbox,
		voxel_size=10.0,
		outmap_fname=os.path.join(data_path, "torus_yeast_fitted.mrc"),
	)

	return


def score_torus_map_with_experimental_map():
	"""
	Score a generated torus map against an experimental reference map.

	This uses the library's scoring function to quantify overlap/similarity between
	two maps (e.g., experimental vs generated).

	Returns:
		None

	Example:
		.. code-block:: python

			fitter = bf.BagelFitter()
			mapfile1 = os.path.join(data_path, "Yeast_C8_Double_MR_center.mrc")
			mapfile2 = os.path.join(data_path, "torus_yeast_fitted.mrc")
			fitter.score_torus_maps(mapfile1, mapfile2)
	"""
	import bagelfit as bf

	# -----------------------------------------------------#
	# Score a torus map against an experimental map
	# -----------------------------------------------------#
	data_path = os.path.join(script_dir, "./yeast_membrane/")
	fitter = bf.BagelFitter()

	mapfile1 = os.path.join(data_path, "Yeast_C8_Double_MR_center.mrc")
	mapfile2 = os.path.join(data_path, "torus_yeast_fitted.mrc")
	fitter.score_torus_maps(mapfile1, mapfile2)

	return


def generate_bestfit_torus_map():
	"""
	Grid-search torus parameters against an experimental map and save the best-fit map.

	This example loads an experimental map, searches over (R, r, theta) ranges, then
	writes the highest-scoring fitted torus map to disk.

	Returns:
		None

	Example:
		.. code-block:: python

			fitter = bf.BagelFitter()
			fitter.load_exprimental_map(os.path.join(data_path, "Yeast_C8_Double_MR_center.mrc"))

			tor_R_range = (660, 670, 10)
			tor_r_range = (140, 160, 20)
			tor_th_range = (55, 65, 10)
			extension = 0.0

			best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)
			fitter.write_torusmap_to_file(os.path.join(data_path, "torus_yeast_fitted.mrc"))
	"""
	import bagelfit as bf

	# -----------------------------------------------------#
	# Fit torus parameters to an experimental nuclear membrane map
	# and write the best-scoring torus map to disk.
	# -----------------------------------------------------#
	data_path = os.path.join(script_dir, "./yeast_membrane/")
	fitter = bf.BagelFitter()

	# Load the experimental map used as the fitting target.
	fitter.load_exprimental_map(os.path.join(data_path, "Yeast_C8_Double_MR_center.mrc"))

	# Parameter ranges for grid search: (start, stop, step) or library-defined convention.
	tor_R_range = (660, 670, 10)
	tor_r_range = (140, 160, 20)
	tor_th_range = (55, 65, 10)

	extension = 0.0

	# Fit a binary torus model (occupancy-based). Use fit_nonbinary_torus for density-based fitting if desired.
	best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)
	# best_torus = fitter.fit_nonbinary_torus(tor_R_range, tor_r_range, tor_th_range, extension)

	fitter.write_torusmap_to_file(os.path.join(data_path, "torus_yeast_fitted.mrc"))

	return


if __name__ == '__main__':

	# ---------------------------------
	# Examples: generate torus maps and/or score them against experimental data
	# ---------------------------------

	# Generate a single binary torus map in a cubic bounding box
	# generate_binary_torus_map(boundingbox=2500, extension=1000)

	# Generate a non-binary torus map using an experimental map context.
	# Note: bounding box dimensions may need to be compatible with the experimental map.
	# generate_nonbinary_torus_map(
	# 	experimental_map='./yeast_membrane/Yeast_C8_Double_MR_center_resampled.mrc',
	# 	boundingbox=3000,
	# 	extension=0,
	# )

	# Score a generated torus map against an experimental map
	# score_torus_map_with_experimental_map()

	# Fit torus parameters and write the best-fit torus map
	# generate_bestfit_torus_map()

	# ---------------------------------
	# Examples: extended membrane around one or more toroids
	# ---------------------------------

	# Example: extension chosen to cover a square domain of side 2500 (diagonal-based estimate)
	# generate_binary_torus_map(boundingbox=2500, extension=1768)  # ~(2*(1250**2))**0.5

	# Generate two toroids in one map
	# generate_two_binary_torus_map()

	# Generate four toroids in one map
	generate_four_binary_torus_map()