import numpy as np
from numpy import linalg as LA
from shapely.geometry import Polygon

class SurfaceGML:
    def __init__(self, gml_surface, boundary=None):
        self.gml_surface = gml_surface
        self.name = boundary
        self.surface_area = None
        self.surface_orientation = None
        self.surface_tilt = None

        try:
            self.surface_area = self.calculate_area()
            self.surface_orientation = self.get_gml_orientation()
            self.surface_tilt = self.get_gml_tilt()
            print(f"area: {self.surface_area}, orientation: {self.surface_orientation}, tilt: {self.surface_tilt}")
        except Exception as e:
            print(f"Error processing surface '{self.name}': {e}")

    def fit_plane(self, points):
        centroid = points.mean(axis=0)
        centered_points = points - centroid
        u, s, vh = LA.svd(centered_points)
        normal = vh[2, :]
        
        if LA.norm(normal) == 0:
            raise ValueError("Degenerate normal vector.")
        
        normal /= LA.norm(normal)
        
        return normal, centroid

    def project_to_2D(self, points, normal, centroid):
        try:
            if np.allclose(normal, [0, 0, 1]):
                in_plane_x = np.array([1, 0, 0])
            else:
                in_plane_x = np.cross(normal, [0, 0, 1])
                if LA.norm(in_plane_x) == 0:
                    in_plane_x = np.array([1, 0, 0])  # Fallback
                else:
                    in_plane_x /= LA.norm(in_plane_x)
            in_plane_y = np.cross(normal, in_plane_x)
            in_plane_y /= LA.norm(in_plane_y)

            translated = points - centroid
            x_coords = np.dot(translated, in_plane_x)
            y_coords = np.dot(translated, in_plane_y)

            return np.vstack((x_coords, y_coords)).T
        except Exception as e:
            print(f"Error projecting points to plane: {e}")
            return None

    def calculate_area(self):
        points = np.array(self.gml_surface).reshape(-1, 3)
        if len(points) < 3:
            print(f"Surface '{self.name}' has fewer than 3 points. Area set to 0.0.")
            return 0.0

        try:
            normal, centroid = self.fit_plane(points)
            projected_2d = self.project_to_2D(points, normal, centroid)

            if projected_2d is None:
                print(f"Projection failed for surface '{self.name}'. Area set to 0.0.")
                return 0.0

            polygon = Polygon(projected_2d)
            if not polygon.is_valid:
                polygon = polygon.buffer(0)
                if not polygon.is_valid:
                    print(f"Surface '{self.name}' has an invalid polygon. Area set to 0.0.")
                    return 0.0

            area = polygon.area
            if area is None:
                area = 0.0
            return area
        except Exception as e:
            print(f"Error calculating area for surface '{self.name}': {e}")
            return 0.0

    def calculate_orientation(self):
        points = np.array(self.gml_surface).reshape(-1, 3)
        normal, _ = self.fit_plane(points)
        

        if np.isclose(normal[2], 1.0, atol=1e-3):
            orientation = -1  # Roof
            return orientation
        elif np.isclose(normal[2], -1.0, atol=1e-3):
            orientation = -2  # Floor
            return orientation
        else:
            phi = np.degrees(np.arctan2(normal[1], normal[0]))
            if phi < 0:
                phi += 360
            teaser_orientation = (90 - phi) % 360
            if teaser_orientation > 180:
                teaser_orientation -= 360
            return teaser_orientation

    def calculate_tilt(self):
        points = np.array(self.gml_surface).reshape(-1, 3)
        normal, _ = self.fit_plane(points)

        tilt = np.degrees(np.arccos(normal[2]))
        
        return tilt

    def _debug(self, message):
        print(message)

    def get_gml_tilt(self):
        """calc the tilt of a gml_surface defined by 4 or 5 gml coordinates

        Surface needs to be planar

        Returns
        ----------
        surface_tilt : float
            returns the orientation of the surface
        """

        gml_surface = np.array(self.gml_surface)
        gml1 = gml_surface[0:3]
        gml2 = gml_surface[3:6]
        gml3 = gml_surface[6:9]

        vektor_1 = gml2 - gml1
        vektor_2 = gml3 - gml1

        normal_1 = np.cross(vektor_1, vektor_2)
        z_axis = np.array([0, 0, 1])

        self.surface_tilt = np.arccos(np.dot(normal_1, z_axis) / (LA.norm(
            z_axis) * LA.norm(normal_1))) * 360 / (2 * np.pi)

        if self.surface_tilt == 180:
            self.surface_tilt = 0.0
        elif str(self.surface_tilt) == "nan":
            self.surface_tilt = 90
        return self.surface_tilt

    def get_gml_orientation(self):
        """calc the orientation of a gml_surface defined by 4 or 5 gml
        coordinates

        Surface needs to be planar, the orientation returned is in TEASER
        coordinates

        Returns
        ----------
        surface_orientation : float
            returns the orientation of the surface
        """

        gml_surface = np.array(self.gml_surface)
        gml1 = gml_surface[0:3]
        gml2 = gml_surface[3:6]
        gml3 = gml_surface[6:9]
        gml4 = gml_surface[9:12]
        if len(gml_surface) > 12:
            vektor_1 = gml2 - gml1
            vektor_2 = gml4 - gml1
        else:
            vektor_1 = gml2 - gml1
            vektor_2 = gml3 - gml1

        normal_1 = np.cross(vektor_1, vektor_2)
        normal_uni = normal_1 / LA.norm(normal_1)
        phi = None
        if normal_uni[0] > 0:
            phi = np.arctan(normal_uni[1] / normal_uni[0])
        elif normal_uni[0] < 0 <= normal_uni[1]:
            phi = np.arctan(normal_uni[1] / normal_uni[0]) + np.pi
        elif normal_uni[0] < 0 > normal_uni[1]:
            phi = np.arctan(normal_uni[1] / normal_uni[0]) - np.pi
        elif normal_uni[0] == 0 < normal_uni[1]:
            phi = np.pi / 2
        elif normal_uni[0] == 0 > normal_uni[1]:
            phi = -np.pi / 2

        if phi is None:
            pass
        elif phi < 0:
            self.surface_orientation = (phi + 2 * np.pi) * 360 / (2 * np.pi)
        else:
            self.surface_orientation = phi * 360 / (2 * np.pi)

        if self.surface_orientation is None:
            pass
        elif 0 <= self.surface_orientation <= 90:
            self.surface_orientation = 90 - self.surface_orientation
        else:
            self.surface_orientation = 450 - self.surface_orientation

        if normal_uni[2] == -1:
            self.surface_orientation = -2
        elif normal_uni[2] == 1:
            self.surface_orientation = -1
        return self.surface_orientation

    def unit_normal(self, a, b, c):
        """calculates the unit normal vector of a surface described by 3 points

        Parameters
        ----------

        a : float
            point 1
        b : float
            point 2
        c : float
            point 3

        Returns
        ----------

        unit_normal : list
            unit normal vector as a list

        """

        x = np.linalg.det([[1, a[1], a[2]],
                           [1, b[1], b[2]],
                           [1, c[1], c[2]]])
        y = np.linalg.det([[a[0], 1, a[2]],
                           [b[0], 1, b[2]],
                           [c[0], 1, c[2]]])
        z = np.linalg.det([[a[0], a[1], 1],
                           [b[0], b[1], 1],
                           [c[0], c[1], 1]])
        magnitude = (x**2 + y**2 + z**2)**.5
        return x / magnitude, y / magnitude, z / magnitude
