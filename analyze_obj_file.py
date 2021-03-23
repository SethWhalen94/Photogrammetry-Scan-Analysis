import pywavefront
from pywavefront import visualization
from sympy import *
import numpy as np
from create_sphere import create_sphere
from find_plane_intersection import find_ray_unit_vector, equation_line_3D
import trimesh
import math
# v_dict = {}
# vn_dict = {}
# vt_dict = {}
# f_dict = {}
# g_dict = {}

# Function to normalize a vector
def normalize(vector):
	return vector / np.linalg.norm(vector)

# find the equation for a plane given 3 vertices (points)
def equation_plane(x1, y1, z1, x2, y2, z2, x3, y3, z3):
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
	print("equation of plane is ",)
	print(a, "x +",)
	print(b, "y +",)
	print(c, "z +",)
	print(d, "= 0.")


# ===========================================
# Opens OBJ file and creates dicts for each component
# RETURNS: a tuple of dictionaries containing OBJ file components
# ===========================================
def read_obj_file(filename):
	v_dict = {}
	vn_dict = {}
	vt_dict = {}
	f_dict = {}
	g_dict = {}
	v = 1
	vn = 1
	vt = 1
	f = 1
	g = 1

	with open(filename, mode='r') as file:
		lines = file.readlines()
		for line in lines:
			# Add vertex to vertex dict
			line = line.strip('\n')
			if line.startswith('v '):
				v_xyz = line.split('v ')[1]
				v_xyz = v_xyz.split(' ')				# xyz is now in for [str(x) str(y-coord) str(z-coord)
				v_xyz = [float(i) for i in v_xyz if i !='']		# now in the form [float(x-coord) float(y-coord) float(z-coord)]
				v_dict[v] = v_xyz
				v += 1
			# Add vertex normals to vertex normals dict
			elif line.startswith('vn '):
				vn_xyz = line.split('vn ')[1]
				vn_xyz = vn_xyz.split(' ')  			# xyz is now in for [str(x) str(y-coord) str(z-coord)
				vn_xyz = [float(i) for i in vn_xyz]  	# now in the form [float(x-coord) float(y-coord) float(z-coord)]
				vn_dict[vn] = vn_xyz
				vn += 1
			# Add vertex textures to vertex textures dict
			elif line.startswith('vt '):
				vt_dict[vt] = line.split('vt ')[1]
				vt += 1
			# Add face to faces dict
			elif line.startswith('f '):
				f_verts = line.split('f ')[1]
				f_verts = f_verts.split(' ')							# f_verts in form [str(v1/vn1/vt1) str(v2/vn2/vt3) ... str(vN/vnN/vtN)]

				if '/' in f_verts[0]:
					f_verts = [int(v.split('/')[0]) for v in f_verts if v != '' ]	# f_verts in form [int(v1) int(v2) ... int(vN)]
				else:
					f_verts = [int(v) for v in f_verts]					# f_verts in form [int(v1) int(v2) ... int(vN)]

				f_dict[f] = f_verts
				f += 1
			# Add group to groups dict
			elif line.startswith('g '):
				g_dict[g] = line.split('g ')[1]
				g += 1

	return v_dict, f_dict, vn_dict, vt_dict


def cartesian_to_Spherical(xyz):
	#takes list xyz (single coord)
	print("coordinate = ", xyz)
	xyz = xyz.split(' ')
	print("coordinate = ", xyz)
	x       = xyz[0]
	y       = xyz[1]
	z       = xyz[2]
	r       =  sqrt(x*x + y*y + z*z)
	theta   =  acos(z/r)*180/ pi #to degrees
	phi     =  atan2(y,x)*180/ pi
	return [r,theta,phi]

# ===========================================
# Converts cartesian to polar radius
# ===========================================
def xyz_to_radius(xyz):
	# takes list xyz (single coord)
	x = xyz[0]
	y = xyz[1]
	z = xyz[2]
	r = sqrt(x * x + y * y + z * z)
	return r

# ===========================================
# Finds the min and max vertices in dictionary
# ===========================================
def get_max_coordinate(v_dict):
	min = max = min_point = max_point = None
	# Find largest radius for spherical coordinates
	for index, vertex in v_dict.items():

		if max == None and min == None and min_point == None and max_point == None:
			max = min = np.linalg.norm(vertex)
			max_point = min_point = index
		elif np.linalg.norm(vertex) > max:
			max = np.linalg.norm(vertex)  # new max value
			max_point = index
		elif np.linalg.norm(vertex) < min:
			min = np.linalg.norm(vertex)  # new min value
			min_point = index

	print("min is ", min)
	print("max is ", max)
	#print("v_dict[min] = ", v_dict[min_point])
	min_point = Point(v_dict[min_point][0], v_dict[min_point][1], v_dict[min_point][2])
	max_point = Point(v_dict[max_point][0], v_dict[max_point][1], v_dict[max_point][2])

	return max, min, max_point, min_point

# ===========================================
# Finds the midpoint of 2 points on a 3D line
# ===========================================
def find_line_midpoint(initial_pt, final_pt):
	slope_x = final_pt.x - initial_pt.x
	slope_y = final_pt.y - initial_pt.y
	slope_z = final_pt.z - initial_pt.z

	slope = np.array([slope_x, slope_y, slope_z])
	position_0 = np.array([initial_pt.x, initial_pt.y, initial_pt.z])

	#equation for 3D line is <x, y, z> = poition_0 + t*slope

	midpoint = position_0 + 0.5 * slope

	midpoint = Point(midpoint[0], midpoint[1], midpoint[2])	# Convert vector to a Point object

	return midpoint

if __name__ == "__main__":

	print("Starting Program...")

	file_name = "Box/model1.obj"

	v_dict, f_dict, vn_dict, vt_dict = read_obj_file(filename=file_name)

	print("Finished gathering OBJ file data...")

	print ("Finding min and max coordinates")
	max, min, max_point, min_point = get_max_coordinate(v_dict)		# get max and min values of vertices

	print(f'The min radius is {min}')
	print(f'The max radius is {max}')

	midpoint = find_line_midpoint(max_point, min_point)		# Find midpoint of max and min points

	ray_origins, ray_directions, num_vertices = create_sphere(file_name="new_sphere1.obj", vertical_lines=20, radius=max, ray_center=midpoint, center_point=midpoint)	# Create a sphere to enclose the 3D object

	# convert to NumPy arrays for trimesh library to use
	ray_origins = np.asarray(ray_origins)
	ray_directions = np.asarray(ray_directions)
	#test_origins = np.array([ray_origins[0], ray_origins[1]])
	#test_directions = np.array([ray_directions[0], ray_directions[1]])

	# start 3D mesh analysis
	mesh = trimesh.load(file_obj=file_name)

	print("Starting ray casting process...")

	locations, index_ray, index_tri = mesh.ray.intersects_location(
		ray_origins=ray_origins,
		ray_directions=ray_directions, multiple_hits=False)

	print("Finished ray casting process...")
	#for i in range(2):
	#	print("Location: ",locations[i])
	#	print("face index: ", index_tri[i])

	# stack rays into line segments for visualization as Path3D
	ray_visualize = trimesh.load_path(np.hstack((
		ray_origins,
		ray_origins + ray_directions)).reshape(-1, 2, 3))

	# =============================================================================
	# Add functionality to find face normals and compare to ray normals
	# Check if locations/index_tri arraya are zero before looking at face normals
	# ==============================================================================

	# unmerge so viewer doesn't smooth
	mesh.unmerge_vertices()
	# make mesh transparent- ish
	#mesh.visual.face_colors = [255, 255, 255, 255]
	#mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
	mesh.visual.face_colors = [100, 100, 100, 100]

	# create a visualization scene with rays, hits, and mesh
	scene = trimesh.Scene([
		mesh,
		ray_visualize,
		trimesh.points.PointCloud(locations)])

	# display the scene

	scene.show()


