import trimesh
import numpy as np
from create_sphere import create_sphere

if __name__ == '__main__':

	# test on a sphere mesh
	mesh = trimesh.load("Folding_Chair/Folding_Chair.obj")

	# create some rays
	ray_origins = np.asarray([[5, 0, 0], [0, 0, -5]])
	ray_directions = np.array([[-1, 0, 0], [0, 0, 1]])

	# Max = 1.904789530958473 radius
	max = 19

	print("finding ray directions")
	#ray_origins, ray_directions, num_vertices = create_sphere(file_name="new_sphere1.obj", vertical_lines=20,
		#													  radius=max)  # Create a sphere to enclose the 3D object

	print("finished finding ray directions")
	"""
	Signature: mesh.ray.intersects_location(ray_origins,
                                            ray_directions,
                                            multiple_hits=True)
    Docstring:
    Return the location of where a ray hits a surface.
    Parameters
    ----------
    ray_origins:    (n,3) float, origins of rays
    ray_directions: (n,3) float, direction (vector) of rays
    Returns
    ---------
    locations: (n) sequence of (m,3) intersection points
    index_ray: (n,) int, list of ray index
    index_tri: (n,) int, list of triangle (face) indexes
    """




	print("starting ray casting")
	# run the mesh- ray test
	locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions, multiple_hits=False)

	print("Finished ray casting")
	# stack rays into line segments for visualization as Path3D
	ray_visualize = trimesh.load_path(np.hstack((
    	ray_origins,
    	ray_origins + ray_directions)).reshape(-1, 2, 3))

	# unmerge so viewer doesn't smooth
	mesh.unmerge_vertices()
	# make mesh transparent- ish
	mesh.visual.face_colors = [255, 255, 255, 255]
	#mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
	# create a visualization scene with rays, hits, and mesh
	scene = trimesh.Scene([
    	mesh,
    	ray_visualize,
    	trimesh.points.PointCloud(locations)])

	# display the scene
	print(locations)
	print(index_tri)
	print(index_ray)
	scene.show()