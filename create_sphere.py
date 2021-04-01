import math
from find_plane_intersection import find_ray_direction, find_ray_unit_vector

class Point():
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	# Member function to return the inverse of point
	def invert_point(self):
		return Point(-self.x, -self.y, -self.z)

DEGS_TO_RAD = 3.14159/180.0
numVertices = 0    # Tallies the number of vertex points added.

#------------------------
#-- Prints a sphere as a "standard sphere" triangular mesh with the specified
#-- number of latitude (nLatitude) and longitude (nLongitude) lines and
#-- writes results to the specified output file (fout).

# PARAMS: Point pt, float radius, int nLatitude, int nLongitude, FILE fout
def printStandardSphere(pt, radius, nLatitude, nLongitude, fout, final_point):

	#global numVertices	# make sure we are referring to GLOBAL numVertices and not a local copy
	numVertices = 0
	vertex_dict = {}	# a vertex dictionary to hold all vertices
	vertex_list = []
	ray_directions_dict = {}
	ray_directions_list = []
	p = s = i = j = 0
	x = y = z = out = float(0)
	nPitch = nLongitude + 1

	pitchInc = (180. / float(nPitch)) * DEGS_TO_RAD
	rotInc   = (360. / float(nLatitude)) * DEGS_TO_RAD

	### PRINT VERTICES:
	# Write vertices to OBJ file !!!!!!
	fout.write("v %g %g %g\n" %(pt.x, pt.y+radius, pt.z) )    # Top vertex.
	numVertices += 1

	ray_origin = Point(pt.x, pt.y + radius, pt.z)
	ray_direction = find_ray_direction(ray_origin, final_pt=final_point)	# Find ray direction
	ray_directions_list.append(ray_direction)									# Add ray direction to list
	ray_directions_dict[numVertices] = ray_direction							# Add ray direction to dictionary
	vertex_dict[numVertices] = [pt.x, pt.y+radius, pt.z]
	vertex_list.append([pt.x, pt.y + radius, pt.z])

	fout.write("v %g %g %g\n" %(pt.x, pt.y-radius, pt.z) )    # Bottom vertex.
	numVertices += 1

	ray_origin = Point(pt.x, pt.y - radius, pt.z)
	ray_direction = find_ray_direction(ray_origin, final_pt=final_point)  	# Find ray direction
	ray_directions_list.append(ray_direction)  									# Add ray direction to list
	ray_directions_dict[numVertices] = ray_direction  							# Add ray direction to dictionary
	vertex_dict[numVertices] = [pt.x, pt.y - radius, pt.z]
	vertex_list.append([pt.x, pt.y - radius, pt.z])

	fVert = numVertices    # Record the first vertex index for intermediate vertices.
	# Generate all "intermediate vertices":
	for p in range(1, nPitch):
		out = radius * math.sin(float(p) * pitchInc)
		if out < 0:
			out = -out    # abs() command won't work with all compilers

		y   = radius * math.cos(p * pitchInc)

		#print("OUT = %g\n", out)    # bottom vertex
		#print("nPitch = %d\n", nPitch)    # bottom vertex
		#for(s=0 s<nLatitude s++)
		for s in range(0, nLatitude):
			x = out * math.cos(s * rotInc)
			z = out * math.sin(s * rotInc)

			fout.write("v %g %g %g\n" %(x + pt.x, y + pt.y, z + pt.z) )
			numVertices += 1

			ray_origin = Point(x + pt.x, y + pt.y, z + pt.z)
			ray_direction = find_ray_direction(ray_origin, final_pt=final_point)  	# Find ray direction
			ray_directions_list.append(ray_direction)  									# Add ray direction to list
			ray_directions_dict[numVertices] = ray_direction  							# Add ray direction to dictionary
			vertex_dict[numVertices] = [x + pt.x, y + pt.y, z + pt.z]
			vertex_list.append([x + pt.x, y + pt.y, z + pt.z])

	### PRINT SQUARE FACES BETWEEN INTERMEDIATE POINTS:

	for p in range(1, nPitch):
		for s in range(0, nLatitude):
			i = p*nLatitude + s
			j = i-nLatitude if s == (nLatitude-1) else i

			# Print face to OBJ file !!!!
			fout.write("f %d %d %d %d\n" %((i+1-nLatitude)+fVert, (j+2-nLatitude)+fVert, (j+2)+fVert, (i+1)+fVert) )

	### PRINT TRIANGLE FACES CONNECTING TO TOP AND BOTTOM VERTEX:

	offLastVerts  = fVert + (nLatitude * (nLongitude-1))

	for s in range(0, nLatitude):
		j = -1 if s == (nLatitude-1) else s

		# Write faces to OBJ file !!!!!!
		fout.write("f %d %d %d\n" %(fVert-1, (j+2)+fVert, (s+1)+fVert) )
		fout.write("f %d %d %d\n" %(fVert, (s+1)+offLastVerts, (j+2)+offLastVerts) )

	return vertex_list, ray_directions_list,  numVertices


# =============================================================================================
# PARAMS: 	file_name = 'some_filename.obj' : STRING
# 		  	Vertical_lines = number of vertical lines in sphere : INT
#			radius = desired radius of the sphere : INT
# RETURNS:	numVertices = Number of vertices created in the sphere OBJ file, vertex_dict
#			vertex_dict = dictionary of vertices from sphere

# This function creates an OBJ file of a sphere given a filename, vertical lines and radius
# =============================================================================================
def create_sphere(file_name, vertical_lines, radius, ray_center = Point(0, 0, 0), center_point = Point(0, 0, 0)):
	nLatitude = vertical_lines  # Number vertical lines.
	nLongitude = int(nLatitude / 2)  # Number horizontal lines.

	# NOTE: for a good sphere use ~half the number of longitude lines than latitude.
	#centerPt = Point(0, 0, 0)  # Position the center of out sphere at (0,0,0).

	fout = open(file_name, "a+")  # open file to append 'a'
	if fout == None:
		raise Exception("Could not open the file specified")

	fout.truncate(0)  # clear contents of file
	vertex_list, ray_directions_list, numVertices = printStandardSphere(pt=center_point, radius=radius, nLatitude=nLatitude, nLongitude=nLongitude, fout=fout, final_point=ray_center)  # Print sphere with radius 10 into file.

	fout.close()  # Close file
	print(f"  # vertices:   {numVertices}\n")
	return vertex_list, ray_directions_list, numVertices


#------------------------
#-- Entry point. This main() function demonstrates how you can
#-- use "printStandardSphere()", but you probably won't
#-- want/need to copy it in your own code.

if __name__ == "__main__":

	create_sphere("python_sphere10.obj", vertical_lines=20, radius=25)





	nLatitude  = 20                  # Number vertical lines.
	nLongitude = int(nLatitude / 2)      # Number horizontal lines.
	# NOTE: for a good sphere use ~half the number of longitude lines than latitude.
	centerPt = Point(0, 0, 0)           # Position the center of out sphere at (0,0,0).

	#if (argc < 2) {
	#  print(stderr, "Must enter: './programname outputfile.obj'\n")
	#  return (-1)
	#}

	#FILE *fout = fopen(argv[1] , "w")
	fout = open("python_sphere11.obj" , "a+")	# open file to append 'a'
	if fout == None:
		raise Exception("Could not open the file specified")

	fout.truncate(0)	# clear contents of file
	numVertices = printStandardSphere(pt=centerPt, radius=25, nLatitude=nLatitude, nLongitude=nLongitude, fout=fout)      # Print sphere with radius 10 into file.

	fout.close()	# Close file
	print(f"  # vertices:   {numVertices}\n")

