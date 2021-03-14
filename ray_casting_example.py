import PIL.Image
import trimesh
import numpy as np
from typing import List


class RayCaster(object):
    def __init__(self, resolution: List, fov: float, mesh: trimesh.Trimesh):
        self.resolution = resolution
        self.fov = fov
        self.mesh = mesh
        self.scene = None
        self._init_scene()

    def _init_scene(self):
        self.scene = self.mesh.scene()
        # any of the automatically generated values can be overridden
        # set resolution, in pixels
        self.scene.camera.resolution = self.resolution
        # set field of view, in degrees
        # make it relative to resolution so pixels per degree is same
        self.scene.camera.fov = self.fov * (self.scene.camera.resolution /
                                            self.scene.camera.resolution.max())

    def cast_ray(self, cam_rot, cam_t):
        """
        :param cam_rot: (3,) array, camera rotation in Euler angle format.
        :param cam_t: (3,) array, camera translation.
        """
        # self.scene.camera.lookat(angles=cam_rot)
        self.scene.set_camera(angles=cam_rot)
        # self.scene.set_camera(angles=cam_rot, distance=cam_t)

        # convert the camera to rays with one ray per pixel
        origin, vectors, pixels = self.scene.camera_rays()

        # intersects_location requires origins to be the same shape as vectors
        origins = np.tile(np.expand_dims(origin, 0), (len(vectors), 1))

        # do the actual ray- mesh queries
        points, index_ray, index_tri = self.mesh.ray.intersects_location(
            origins, vectors, multiple_hits=False)

        # for each hit, find the distance along its vector
        # you could also do this against the single camera Z vector
        depth = trimesh.util.diagonal_dot(points - origin, vectors[index_ray])

        # find the angular resolution, in pixels per radian
        ppr = self.scene.camera.resolution / np.radians(self.scene.camera.fov)
        # convert rays to pixel locations
        angles = self.scene.camera.angles()
        pixel = (angles * ppr).round().astype(np.int64)
        # make sure we are in the first quadrant
        pixel -= pixel.min(axis=0)
        # find pixel locations of actual hits
        pixel_ray = pixel[index_ray] - 1

        # create a numpy array we can turn into an image
        # doing it with uint8 creates an `L` mode greyscale image
        a = np.zeros(self.scene.camera.resolution, dtype=np.uint8)

        # scale depth against range (0.0 - 1.0)
        depth_float = ((depth - depth.min()) / depth.ptp())

        # convert depth into 0 - 255 uint8
        depth_int = (depth_float * 255).astype(np.uint8)
        # assign depth to correct pixel locations
        a[pixel_ray[:, 0], pixel_ray[:, 1]] = depth_int
        return a


if __name__ == '__main__':
    # Conclusion here:
    # There are strange black stides(zero value) within rendered depth image, either with pyembree or not.
    # No similar issue was found.
    from matplotlib import pyplot as plt
    # test on a simple mesh
    m = trimesh.load('compound_bow.ply')
    # m = trimesh.load('/mnt/sdj/dataset/shapenet/shapenet-core/ShapeNetCore.v2/03001627/'
    #                  'eb3029393f6e60713ae92e362c52d19d/models/model_normalized.obj')

    # m1 = next(iter(m.geometry.values()))
    # m1 = next(iter(m.geometry.values()))

    rc = RayCaster([480, 640], 60, m)
    cam_rot = np.array([0, 0, 0])
    cam_t = np.array([0, 0, 0])
    depth = rc.cast_ray(cam_rot, cam_t)
    plt.imshow(depth, cmap='jet')
    plt.show()