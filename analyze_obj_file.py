import pywavefront
from pywavefront import visualization
from sympy import *
import numpy as np
from create_sphere import create_sphere
from find_plane_intersection import find_ray_unit_vector, equation_line_3D
import trimesh
import math
import os

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

# =======================================================================================
# Finds unit vector normal of a given face plane
# PARAM: Face is an entry from the face dictionary, it is an array of vertex indices
# =============================================================================================
def find_face_normal(face, vertex_dict):
	normal = [0, 0, 0]	# initialize vectors


	p1 = np.array(vertex_dict[face[0]])		# Point 1 on face
	p2 = np.array(vertex_dict[face[1]])		# Point 2 on face
	p3 = np.array(vertex_dict[face[2]])		# Point 3 on face

	v1_2 = p2 - p1							# vector from p1 to p2
	v2_3 = p3 - p2							# vector from p2 to p3
	v3_1 = p1 - p3							# vector from p3 to p1

	norm1 = np.cross(v1_2, v2_3)			# calculate normals
	norm2 = np.cross(v2_3, v3_1)
	norm3 = np.cross(v3_1, v1_2)

	# Normalize vectors
	norm1 = normalize(norm1)
	norm2 = normalize(norm2)
	norm3 = normalize(norm3)

	face_normal = (norm1 + norm2 + norm3)/3		# calculate overall face normal

	return face_normal

# ===================================================
# method to find missing image locations
# PARAMS: ray_origins, ray_directions, face_indices, face_dict, vertex_dict
# RETURNS: List of [x, y, z] coordinates
# ===================================================
def find_missing_image_coordiantes(index_ray, ray_origins, ray_directions, face_indices, face_dict, vertex_dict, locations):

	coords_list = []
	ray_origins_list = []
	ray_directions_list = []
	locations_list = []
	# Use index_ray to determine which ray is paired with the current face
	for index,face in enumerate(face_indices):
		face_normal = find_face_normal(face_dict[face], vertex_dict)

		#dot_prod = np.dot(ray_directions[index], face_normal)
		dot_prod = np.dot(ray_directions[index_ray[index]], face_normal)

		if dot_prod >= 0:
															# Normals are pointing in the same direction
															# this means we are missing image data, since we are hitting the 'back' of the face
			coords_list.append(ray_origins[index_ray[index]])			# Add ray_origin to coordinates list
			ray_origins_list.append(ray_origins[index_ray[index]].tolist())
			ray_directions_list.append(ray_directions[index_ray[index]].tolist())
			locations_list.append(locations[index].tolist())

		elif dot_prod < 0:
			# this means we are hitting the face head on, which is correct
			pass
		else:
			raise Exception("There is a problem with the ray/face dot product calculation")


	ray_origins_list, ray_directions_list = compare_missing_image_coordinates(ray_origins_list, ray_directions_list)			# check if we have already seen these ray hits

	# Convert lists to numpy arrays
	locations_list = np.asarray(locations_list)

	return coords_list, ray_origins_list, ray_directions_list, locations_list


# ==========================================================================================
# Function to check existing missing image coordinates with newly found ones so we do not
# re-use already found coordinates
# PARAMS: image_coordinate - coordinates found from current iteration
# ==========================================================================================
def compare_missing_image_coordinates(image_coordinates, camera_directions):

	new_coordinates = []
	new_directions = []
	exists = False

	with open("missing_image_coordinates.txt", 'a+') as coords:




		coords.seek(0)							# move to beginning of file to read currrent contents
		lines = coords.readlines()
		coords.seek(os.SEEK_END)

		for coord, direction in zip(image_coordinates, camera_directions):

			for line in lines:

				curr_coord = line.split('+')[0]

				if str(coord) == curr_coord:
					exists = True
					break

			if not exists:
				coords.write(str(coord) + '+' + str(direction))  		# Add camera location and direction to the coordinates file
				coords.write('\n')
				new_coordinates.append(coord)  		# add it to new coordinate list
				new_directions.append(direction)	# Add it to new direction list
				exists = False


	return np.asarray(new_coordinates), np.asarray(new_directions)			# return numpy array of new coordinates


# ======================================
# Function to clear missing image coordinates file if this is a new scan
# ==========================================================================
def check_new_scan():
	# Check if user wants to clear current file which holds coordinates
	choice = ''
	while choice.lower() != 'n' and choice.lower() != 'e':
		choice = input(
			"Is this for a new or existing scan? Enter: (N new scan, E for existing scan): ")

	if choice.lower() == 'n':
		with open("missing_image_coordinates.txt", 'a+') as coords:
			coords.truncate(0)

# ============================================================
# MAIN method: This contains the driver code for the program
# ============================================================
def main():
	file_name = "Dice/Dice_hole.obj"
	sphere_obj_name = "sphere1.obj"


	check_new_scan()			# Check if this is a new scan, if so we will clear the coordinates file

	v_dict, f_dict, vn_dict, vt_dict = read_obj_file(filename=file_name)

	print("Finished gathering OBJ file data...")

	print("Finding min and max coordinates")
	max, min, max_point, min_point = get_max_coordinate(v_dict)  # get max and min values of vertices

	print(f'The min radius is {min}')
	print(f'The max radius is {max}')

	midpoint = find_line_midpoint(max_point, min_point)  # Find midpoint of max and min points

	print("creating bounding sphere...")
	ray_origins, ray_directions, num_vertices = create_sphere(file_name=sphere_obj_name, vertical_lines=20, radius=max,
															  ray_center=midpoint,
															  center_point=midpoint)  # Create a sphere to enclose the 3D object

	print("Finished bounding sphere...")
	# convert to NumPy arrays for trimesh library to use
	ray_origins = np.asarray(ray_origins)
	ray_directions = np.asarray(ray_directions)

	# start 3D mesh analysis
	mesh = trimesh.load(file_obj=file_name)

	print("Starting ray casting process...")

	locations, index_ray, index_tri = mesh.ray.intersects_location(  # locations = [x, y, z] of ray hit on mesh
		ray_origins=ray_origins,  # index_ray = array or ray indices for hit locations
		ray_directions=ray_directions, multiple_hits=False)

	print("Finished ray casting process...")

	# =============================================================================
	# find face normals and compare to ray normals
	# Check if locations/index_tri arrays are zero before looking at face normals
	# ==============================================================================
	print("Finding missing image coordinates...")
	missing_image_coords, missing_image_ray_origins, missing_image_ray_directions, hit_locations = find_missing_image_coordiantes(
		index_ray=index_ray,
		ray_origins=ray_origins,
		ray_directions=ray_directions,
		face_indices=index_tri,
		face_dict=f_dict,
		vertex_dict=v_dict,
		locations=locations)

	# stack rays into line segments for visualization as Path3D


	print("The length of missing image coordinates is ", len(missing_image_coords))

	# Use this to plot all rays that were cast and the associated hits
	ray_visualize_all = trimesh.load_path(np.hstack((
		ray_origins,
		ray_origins + ray_directions*0.3)).reshape(-1, 2, 3))

	# unmerge so viewer doesn't smooth
	mesh.unmerge_vertices()
	# make mesh transparent- ish
	# mesh.visual.face_colors = [255, 255, 255, 255]
	# mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
	mesh.visual.face_colors = [100, 100, 100, 100]

	# create a visualization scene with rays, hits, and mesh
	# Scene for all rays
	scene_all = trimesh.Scene([
		mesh,
		ray_visualize_all,
		trimesh.points.PointCloud(locations)])

	# Scene for only missing image locations and associated rays
	if len(hit_locations) > 0 and len(missing_image_ray_origins) > 0:
		# Use this to plot only missing image location rays

		ray_visualize_missing_images = trimesh.load_path(np.hstack((
			missing_image_ray_origins,
			missing_image_ray_origins + missing_image_ray_directions)).reshape(-1, 2, 3))

		scene_missing_images = trimesh.Scene([
			mesh,
			ray_visualize_missing_images,
			trimesh.points.PointCloud(hit_locations)])

		scene_missing_images.show()

	else:
		#display the scene
		scene_all.show()


# ===================================================
# MAIN Method
# ===================================================
if __name__ == "__main__":

	print("Starting Program...")
	main()



