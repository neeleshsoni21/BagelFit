"""Summary
"""
from numpy import sqrt

class Torus:

	"""
	Represents a toroidal shape characterized by its major radius (R), minor radius (r), and thickness.
	
	Attributes:
		R (float): Major radius of the torus.
		r (float): Minor radius of the torus.
		thickness (float): Thickness of the torus.
		extension (float): Optional extension parameter (default is 0.0).
		eps (float): Small constant to avoid division by zero errors.
		dmap (Any): Data map for storing computed values (if applicable).
	"""
	
	def __init__(self, R: float, r: float, thickness: float, extension: float = 0.0, center: tuple[float, float, float] = (0.0, 0.0, 0.0)):
		"""
		Initializes a torus with given parameters.
		
		Args:
			R (float): Major radius of the torus.
			r (float): Minor radius of the torus.
			thickness (float): Thickness of the torus.
			extension (float, optional): Extension parameter (default is 0.0).
		"""
		self.R = R
		self.r = r
		self.thickness = thickness
		self.extension = extension
		self.eps = 1E-9
		self.dmap = None

		self.cx, self.cy, self.cz = center

	@property
	def center(self) -> tuple[float, float, float]:
		return (self.cx, self.cy, self.cz)

	def _to_local(self, x: float, y: float, z: float) -> tuple[float, float, float]:
		"""Convert world coords -> torus-local coords."""
		return x - self.cx, y - self.cy, z - self.cz

	def distance(self, x: float, y: float, z: float, d2_xy: float) -> float:
		"""
		Computes the shortest distance from a point (x, y, z) to the torus.
		
		Args:
			x (float): X-coordinate of the point.
			y (float): Y-coordinate of the point.
			z (float): Z-coordinate of the point.
			d2_xy (float): Squared distance from the torus center in the xy-plane.
		
		Returns:
			float: Shortest distance from the given point to the torus.
		"""
		
		d_xy = sqrt(d2_xy)
		
		if d_xy <= self.R:
			d_tx = x - x/d_xy * self.R if d_xy > self.eps else x - self.R
			d_ty = y - y/d_xy * self.R if d_xy > self.eps else y
			dz = z
		else:
			dz = abs(z) - self.r

		return sqrt(dz**2 + d_tx**2 + d_ty**2)

	def contains_point(self, x: float, y: float, z: float) -> bool:
		"""
		Determines whether a point (x, y, z) lies within the toroidal volume.
		
		Args:
			x (float): X-coordinate of the point.
			y (float): Y-coordinate of the point.
			z (float): Z-coordinate of the point.
		
		Returns:
			bool: True if the point is inside the torus, False otherwise.
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
	

	def contains_point2(self, x: float, y: float, z: float) -> bool:
		"""
		Alternative method to check if a point (x, y, z) lies within the toroidal volume.
		
		Args:
			x (float): X-coordinate of the point.
			y (float): Y-coordinate of the point.
			z (float): Z-coordinate of the point.
		
		Returns:
			bool: True if the point is inside the torus, False otherwise.
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

	


