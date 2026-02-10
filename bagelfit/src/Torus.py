"""Summary
"""
from numpy import sqrt

class Torus:

	"""
	Geometric torus used to model the nuclear membrane.
	
	The membrane is represented as a *shell* within the torus tube: a voxel is considered
	"inside the membrane" when it lies within a band of thickness `thickness` near the
	outer surface of the tube (i.e., within the minor radius `r`).
	
	Terminology:
	- R: major radius (distance from torus center to the tube centerline in the XY plane)
	- r: minor radius (tube radius)
	- thickness: membrane-shell thickness within the tube (0 < thickness <= r)
	- extension: optional in-plane (XY) radial allowance used by map generation to include
	  extended membrane beyond the torus ring (ρ_xy <= R + extension), where ρ_xy is the
	  distance to the torus center in the XY plane.
	"""
	
	def __init__(self, R, r, thickness, extension= 0.0, center= (0.0, 0.0, 0.0)):
		"""
		Initialize a torus.
		
		Args:
		    R: Major radius of the torus.
		    r: Minor radius (tube radius).
		    thickness: Membrane-shell thickness within the tube (radial band thickness).
		    extension: Optional in-plane (XY) radial allowance used by map generation (default 0.0).
		    center: Torus center (cx, cy, cz) in the same coordinate system as map voxels.
		"""
		self.R = R
		self.r = r
		self.thickness = thickness
		self.extension = extension
		self.eps = 1E-9
		self.dmap = None

		self.cx, self.cy, self.cz = center

	@property
	def center(self):
		return (self.cx, self.cy, self.cz)

	def _to_local(self, x, y, z):
		"""Convert world coords -> torus-local coords."""
		return x - self.cx, y - self.cy, z - self.cz

	def distance(self, x, y, z, d2_xy):
		"""
		Compute an approximate shortest distance from a point to the torus surface.
		
		Notes:
		    This implementation assumes the input coordinates are already in *torus-local*
		    coordinates (i.e., after subtracting the torus center). The caller may pass the
		    precomputed squared XY radius (d2_xy) for efficiency.
		
		Args:
		    x, y, z: Torus-local coordinates of the point.
		    d2_xy: Squared radius in the XY plane (x^2 + y^2).
		
		Returns:
		    Distance (float) to the torus surface (implementation-specific approximation).
		"""
		
		d_xy = sqrt(d2_xy)
		
		if d_xy <= self.R:
			d_tx = x - x/d_xy * self.R if d_xy > self.eps else x - self.R
			d_ty = y - y/d_xy * self.R if d_xy > self.eps else y
			dz = z
		else:
			dz = abs(z) - self.r

		return sqrt(dz**2 + d_tx**2 + d_ty**2)

	def contains_point(self, x, y, z):
		"""
		Return 1 if the point lies within the torus *membrane shell*, else 0.
		
		A point is considered inside when:
		- it lies within the torus tube (minor radius constraint), and
		- it falls within the shell band defined by `thickness` near the tube boundary.
		
		Args:
		    x, y, z: World coordinates of the point.
		
		Returns:
		    int: 1 if inside the membrane shell, otherwise 0.
		"""

		# convert to torus-local coordinates
		x, y, z = self._to_local(x, y, z)
	
		if abs(z) > self.r:  # Early exit
			return 0

		d2_xy = x**2 + y**2; d2_z = z**2
		RR_square = d2_xy + d2_z
		RR = sqrt(RR_square)

		if ((self.R - self.r)**2 <= RR_square <= self.R**2):

			return int(self.r - self.thickness < self.distance(x, y, z, d2_xy) <= self.r)
		
		else:
		
			return int(RR > self.R and abs(z) >= self.r - self.thickness)
	

	def contains_point2(self, x, y, z):
		"""
		Alternative membrane-shell membership test (kept for reference).
		
		This method is not used by the main fitting pipeline and may be inconsistent with
		`contains_point()`.
		
		Args:
		    x, y, z: World coordinates of the point.
		
		Returns:
		    int: 1 if inside the membrane shell, otherwise 0.
		"""
		if abs(z) <= self.r:
			RR = sqrt(x**2 + y**2 + z**2)
			if (RR >= self.R - self.r) and (RR <= self.R):
				
				dd = self.distance(x, y, z)
				if dd <= self.r  and dd > self.r - self.thickness:
					return 1
				else:
					return 0
			elif RR > self.R and abs(z)>=self.r - self.thickness:
				return 1
			else:
				return 0
			
		else:
			return 0

	


