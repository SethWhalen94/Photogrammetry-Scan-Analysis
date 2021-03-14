import pywavefront
from pywavefront import visualization
from sympy import *
import numpy as np
import math
v_dict = {}
vn_dict = {}
vt_dict = {}
f_dict = {}
g_dict = {}

# Function to normalize a vector
def normalize(vector):
    return vector / np.linalg.norm(vector)

# Function to determine the intersection of a ray and a sphere
def sphere_intersect(center, radius, ray_origin, ray_direction):
    b = 2 * np.dot(ray_direction, ray_origin - center)
    c = np.linalg.norm(ray_origin - center) ** 2 - radius ** 2
    delta = b ** 2 - 4 * c
    if delta > 0:
        t1 = (-b + np.sqrt(delta)) / 2
        t2 = (-b - np.sqrt(delta)) / 2
        if t1 > 0 and t2 > 0:
            return min(t1, t2)
    return None

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

def read_obj_file(filename):
	v = 1
	vn = 1
	vt = 1
	f = 1
	g = 1

	with open(filename, mode='r') as file:
		lines = file.readlines()
		for line in lines:
			if line.startswith('v '):
				v_xyz = line.split('v ')[1]
				v_xyz = v_xyz.split(' ')				# xyz is now in for [str(x) str(y-coord) str(z-coord)
				v_xyz = [float(i) for i in v_xyz]		# now in the form [float(x-coord) float(y-coord) float(z-coord)]
				v_dict[v] = v_xyz
				v += 1
			elif line.startswith('vn '):
				vn_xyz = line.split('vn ')[1]
				vn_xyz = vn_xyz.split(' ')  			# xyz is now in for [str(x) str(y-coord) str(z-coord)
				vn_xyz = [float(i) for i in vn_xyz]  	# now in the form [float(x-coord) float(y-coord) float(z-coord)]
				vn_dict[vn] = vn_xyz
				vn += 1
			elif line.startswith('vt '):
				vt_dict[vt] = line.split('vt ')[1]
				vt += 1
			elif line.startswith('f '):
				f_dict[f] = line.split('f ')[1]
				f += 1
			elif line.startswith('g '):
				g_dict[g] = line.split('g ')[1]
				g += 1


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

def xyz_to_radius(xyz):
	# takes list xyz (single coord)
	x = xyz[0]
	y = xyz[1]
	z = xyz[2]
	r = sqrt(x * x + y * y + z * z)
	return r

if __name__ == "__main__":
	read_obj_file('Wooden chair.obj')

	# for index in range(10):
	#    print(f'{index+1}: {v_dict[index+1]}\n')
	# for index in range(10):
	#    print(f'{index+1}: {vn_dict[index+1]}\n')
	# for index in range(10):
	#    print(f'{index+1}: {vt_dict[index+1]}\n')
	# for index in range(10):
	#    print(f'{index+1}: {f_dict[index+1]}\n')

	max = 0
	point = 0

	#	Find largest radius for spherical coordinates
	for index, vertex in v_dict.items():
		if xyz_to_radius(vertex)  > max:
			max = xyz_to_radius(vertex)		# new max value
			point = index

	print(f'The greatest radius point is {v_dict[point]}')
	print(f'The radius is {max}')