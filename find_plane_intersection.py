import numpy as np
import math


class Point():
	def __init__(self, x=0.0, y=0.0, z=0.0):
		self.x = x
		self.y = y
		self.z = z

	# Member function to return the inverse of point
	def invert_point(self):
		return Point(-self.x, -self.y, -self.z)

vertex_dict = {}
face_dict = {}
sphere_points = {}
intersections = []
intersection_normals = []

def normalize(vector):
	return vector / np.linalg.norm(vector)


# ================================================
# Method to find coefficients for plane equation for a given face
# PARAM: face = [vert1, vert2, vert3, ...]
# ================================================
def equation_plane(face):

	x1, y1, z1 = vertex_dict[face[0]][0], vertex_dict[face[0]][1], vertex_dict[face[0]][2]	# Gets x, y, z values of first vertex in face
	x2, y2, z2 = vertex_dict[face[1]][0], vertex_dict[face[1]][1], vertex_dict[face[1]][2]  # Gets x, y, z values of second vertex in face
	x3, y3, z3 = vertex_dict[face[2]][0], vertex_dict[face[2]][1], vertex_dict[face[2]][2]  # Gets x, y, z values of third vertex in face

	a1 = x2 - x1
	b1 = y2 - y1
	c1 = z2 - z1
	a2 = x3 - x1
	b2 = y3 - y1
	c2 = z3 - z1
	a = b1 * c2 - b2 * c1
	b = a2 * c1 - a1 * c2
	c = a1 * b2 - b1 * a2
	d = (- a * x1 - b * y1 - c * z1)
	# print("equation of plane is ",)
	# print(a, "x +",)
	# print(b, "y +",)
	# print(c, "z +",)
	# print(d, "= 0.")

	return {'a': a, 'b': b, 'c': c, 'd': d}	# Return dictionary with constants

# ==================================================================================================
#	Equation for 3D line is <x, y, z> = <x0, y0, z0> + t<Mx, My, Mz> where M is slope, t is scalar
# ==================================================================================================
def equation_line_3D(initial_pt, final_pt):
	slope_x = final_pt.x - initial_pt.x
	slope_y = final_pt.y - initial_pt.y
	slope_z = final_pt.z - initial_pt.z

	slope = np.array([slope_x, slope_y, slope_z])
	position_0 = np.array([initial_pt.x, initial_pt.y, initial_pt.z])

	#equation for 3D line is <x, y, z> = poition_0 + t*slope

	return position_0, slope

# ==================================================================================================
#	finds direction of ray and returns the unit vector for said direction
# 	PARAMS: initial_pt, final_pt : type POINT
# ==================================================================================================
def find_ray_direction(initial_pt, final_pt):
	slope_x = final_pt.x - initial_pt.x
	slope_y = final_pt.y - initial_pt.y
	slope_z = final_pt.z - initial_pt.z

	slope = [float(slope_x), float(slope_y), float(slope_z)]
	return slope/np.linalg.norm(slope)


# ===================================
# check if a point lies on a plane
# ===================================
def plane_intersect(face, plane_coeff, point):


	# Get min and max bounds of the face plane
	min_x, max_x, min_y, max_y, min_z, max_z = min_max_of_face(face, vertex_dict)

	# Check if point is outside of the planes bounds
	point_in_plane = point_in_range(point, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y, min_z=min_z, max_z=max_z)

	# Calculate plane equation using x y z of point
	if point_in_plane:
		equation = plane_coeff['a']*point[0] + plane_coeff['b']*point[1] + plane_coeff['c']*point[2] + plane_coeff['d']
		intersects = True if math.isclose(equation, 0, abs_tol=0.00015) else False
		return intersects
	else:
		# Point is not in planes range
		return False


# ==========================================================
# Method to find min and max x, y, z coordinates of a face
# ==========================================================
def min_max_of_face(face, vertices):

	min_x = min_y = min_z = float('inf')
	max_x = max_y = max_z = 0

	for vert in face:

		curr_xyz = vertices[vert]

		# Check x coordinate against current min_x and max_x
		min_x = curr_xyz[0] if curr_xyz[0] < min_x else min_x
		max_x = curr_xyz[0] if curr_xyz[0] > max_x else max_x

		# Check y coordinate against current min_y and max_y
		min_y = curr_xyz[1] if curr_xyz[1] < min_y else min_y
		max_y = curr_xyz[1] if curr_xyz[1] > max_y else max_y

		# Check z coordinate against current min_z and max_z
		min_z = curr_xyz[2] if curr_xyz[2] < min_z else min_z
		max_z = curr_xyz[2] if curr_xyz[2] > max_z else max_z

	return min_x, max_x, min_y, max_y, min_z, max_z

# ==========================================================
# Method to find min and max x, y, z coordinates in vertex dictionary
# ==========================================================
def min_max_of_vertices(vertices):

	min_x = min_y = min_z = float('inf')
	max_x = max_y = max_z = 0

	for vert in vertices:

		curr_xyz = vertices[vert]	# in form [x, y, z]

		# Check x coordinate against current min_x and max_x
		min_x = curr_xyz[0] if curr_xyz[0] < min_x else min_x
		max_x = curr_xyz[0] if curr_xyz[0] > max_x else max_x

		# Check y coordinate against current min_y and max_y
		min_y = curr_xyz[1] if curr_xyz[1] < min_y else min_y
		max_y = curr_xyz[1] if curr_xyz[1] > max_y else max_y

		# Check z coordinate against current min_z and max_z
		min_z = curr_xyz[2] if curr_xyz[2] < min_z else min_z
		max_z = curr_xyz[2] if curr_xyz[2] > max_z else max_z

	return min_x, max_x, min_y, max_y, min_z, max_z

# ==========================================================
# Method to check if a point is in a x, y, z range
# ==========================================================
def point_in_range(point, min_x, max_x, min_y, max_y, min_z, max_z, tolerance=0.00001):

	# Check if point is outside of the planes bounds, or in case of min = max, check if point is within a tolerance of 0.00001
	if not min_x <= point[0] <= max_x and not(math.isclose(point[0], min_x, rel_tol=tolerance) and math.isclose(point[0], max_x, rel_tol=tolerance)):
		return False
	if not min_y <= point[1] <= max_y and not(math.isclose(point[1], min_y, rel_tol=tolerance) and math.isclose(point[1], max_y, rel_tol=tolerance)):
		return False
	if not min_z <= point[2] <= max_z and not(math.isclose(point[2], min_z, rel_tol=tolerance) and math.isclose(point[2], max_z, rel_tol=tolerance)):
		return False

	return True


# =========================================
# Method to calculate point on 3D line
# <x, y, z> = <x0, y0, z0> + t<Mx, My, Mz>
# =========================================
def point_on_3D_line(initial_point, slope, scalar):

	point = initial_point + scalar*slope

	return point

# =================================================
# Returns a value rounded up to a specific number of decimal places.
# =================================================
def round_decimals_down(number:float, decimals:int=5):

	if not isinstance(decimals, int):
		raise TypeError("decimal places must be an integer")
	elif decimals < 0:
		raise ValueError("decimal places has to be 0 or more")
	elif decimals == 0:
		return math.floor(number)

	factor = 10 ** decimals
	return math.floor(number * factor) / factor


# ====================================================================
# Function to check which faces fall relatively inside the rays path
# PARAMS: faces - dictionary, ray_initial - initial Point of ray, ray_final - final Point of ray
# ====================================================================
def parse_faces(faces, ray_initial, ray_final):
	valid_faces = []
	min_x = min(ray_initial.x, ray_final.x)
	max_x = max(ray_initial.x, ray_final.x)
	min_y = min(ray_initial.y, ray_final.y)
	max_y = max(ray_initial.y, ray_final.y)
	min_z = min(ray_initial.z, ray_final.z)
	max_z = max(ray_initial.z, ray_final.z)

	for index, face in faces.items():
		face_min_x, face_max_x, face_min_y, face_max_y, face_min_z, face_max_z = min_max_of_face(face,
																								 vertices=vertex_dict)

		# min_x, min_y, min_z
		if face_min_x in range(min_x, max_x) and face_min_y in range(min_y, max_y) and face_min_z in range(min_z,
																										   max_z):
			valid_faces.append(face)

		# min_x, min_y, max_z
		elif face_min_x in range(min_x, max_x) and face_min_y in range(min_y, max_y) and face_max_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

		# min_x, max_y, min_z
		elif face_min_x in range(min_x, max_x) and face_max_y in range(min_y, max_y) and face_min_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

		# min_x, max_y, max_z
		elif face_min_x in range(min_x, max_x) and face_max_y in range(min_y, max_y) and face_max_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

		# max_x, min_y, min_z
		elif face_max_x in range(min_x, max_x) and face_min_y in range(min_y, max_y) and face_min_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

		# max_x, min_y, max_z
		elif face_max_x in range(min_x, max_x) and face_min_y in range(min_y, max_y) and face_max_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

		# max_x, max_y, min_z
		elif face_max_x in range(min_x, max_x) and face_max_y in range(min_y, max_y) and face_min_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

		# max_x, max_y, max_z
		elif face_max_x in range(min_x, max_x) and face_max_y in range(min_y, max_y) and face_max_z in range(min_z,
																											 max_z):
			valid_faces.append(face)

	return valid_faces

# =================================================
# Casts a ray from initial point to final point
# checks if ray intersects a face plane
# =================================================
def cast_ray(face_dict, initial_point, final_point):
	# Ray casting path
	ray_initial = initial_point
	ray_final = final_point

	# 3D line parameters
	line_initial, line_slope = equation_line_3D(ray_initial, ray_final)

	# steps for ray being cast,
	steps = np.linspace(0, 1, num=1000000)

	vert_min_x, vert_max_x, vert_min_y, vert_max_y, vert_min_z, vert_max_z = min_max_of_vertices(vertex_dict)
	# print(steps)
	for step in steps:


		point = point_on_3D_line(line_initial, line_slope, step)  # Calculate point on 3D line

		# check if point is within min, max values of vertex values
		point_valid = point_in_range(point, min_x=vert_min_x, max_x=vert_max_x, min_y=vert_min_y, max_y=vert_max_y,
									 min_z=vert_min_z, max_z=vert_max_z, tolerance=0.00001)

		# go through all faces in dictionary if point is valid
		if point_valid:
			for face_number in face_dict:
				# Check if point intersects faces plane
				plane_coeff = equation_plane(face_dict[face_number])
				intersects = plane_intersect(face=face_dict[face_number], plane_coeff=plane_coeff, point=point)
				if intersects:
					intersections.append(point)
					intersection_normals.append(normalize(line_slope))


# =================================================
# Finds unit vector normal of a given face plane
# PARAM: Face is an entry from the face dictionary, it is an array of vertex indices
# =================================================
def find_face_normal(face):
	normal = [0, 0, 0]	# initialize vectors

	for i in range(1, len(face)):
		vect_a = np.array(vertex_dict[i])
		for j in range(i + 1, len(face)):
			vect_b = np.array(vertex_dict[j])
			normal = np.cross(vect_a, vect_b)

			if not np.all((normal == 0)):	# inner loop: check if normal in non-zero
				break
		if not np.all((normal == 0)):		# outer loop: check if normal in non-zero
			break

	#print("vect_a is ", vect_a)
	#print("vect_b is ", vect_b)
	#print("Normal vector is ", normal)
	normal = normalize(normal)  # Normalize vector so it is a unit vector

	return normal

# ===========================================
# Take the slope vector of a ray and calculates its unit vector
# ============================================
def find_ray_unit_vector(ray_slope):
	return normalize(ray_slope)

# =================================================
# Creates a test sphere with 14 coordinate points
# =================================================
def create_circle_dict():
	sphere_dict = {}
	sphere_dict[1] = [3.0, 0.0, 0.0]
	sphere_dict[2] = [-3.0, 0.0, 0.0]

	sphere_dict[3] = [0.0, 3.0, 0.0]
	sphere_dict[4] = [0.0, -3.0, 0.0]

	sphere_dict[5] = [0.0, 0.0, 3.0]
	sphere_dict[6] = [0.0, 0.0, -.0]

	sphere_dict[7] = [3.0, 3.0, 3.0]
	sphere_dict[8] = [3.0, 3.0, -3.0]

	sphere_dict[9] = [3.0, -3.0, 3.0]
	sphere_dict[10] = [3.0, -3.0, 3.0]

	sphere_dict[11] = [-3.0, 3.0, 3.0]
	sphere_dict[12] = [-3.0, 3.0, -3.0]

	sphere_dict[13] = [-3.0, -3.0, 3.0]
	sphere_dict[14] = [-3.0, -3.0, 3.0]

	return sphere_dict



if __name__ == "__main__":
	vertex_dict[1] = [2.0324, -3.0, -3.0]
	vertex_dict[2] = [2.0324, 3.0, -3.0]
	vertex_dict[3] = [2.0324, 3.0, 3.0]
	vertex_dict[4] = [2.0324, -3.0, 3.0]

	sphere_dict = create_circle_dict()
	face_dict[1] = [1, 2, 3, 4]
	ray = Point(1.0, 0, 0.0)
	# for each vertex in sphere vertices
	# set initial point = vertex, and final_point = -vertex
	for index,ray in sphere_dict.items():

		ray = Point(ray[0], ray[1], ray[2])	# create point with xyz coordinates of vertex

		cast_ray(face_dict=face_dict, initial_point=ray, final_point=ray.invert_point())

	print(intersection_normals)