import numpy as np
from numpy import linalg as LA
from shapely.geometry import Polygon

class SurfaceGML:
    """
    A class to represent and analyze a geometric surface defined by GML coordinates.

    Attributes
    ----------
    gml_surface : list
        A list of coordinates defining the surface.
    name : str, optional
        The name or boundary identifier of the surface.
    surface_area : float
        The calculated area of the surface.
    surface_orientation : float
        The orientation of the surface in TEASER coordinates.
    surface_tilt : float
        The tilt of the surface in degrees.

    Methods
    -------
    fit_plane(points):
        Fits a plane to the given 3D points and returns the normal vector and centroid.
    
    project_to_2D(points, normal, centroid):
        Projects the 3D points onto a 2D plane defined by the normal vector and centroid.
    
    calculate_area():
        Calculates the area of the surface based on the GML coordinates.
    
    calculate_orientation():
        Calculates the orientation of the surface based on the normal vector.
    
    calculate_tilt():
        Calculates the tilt of the surface based on the normal vector.
    
    get_gml_tilt():
        Calculates the tilt of the GML surface defined by its coordinates.
    
    get_gml_orientation():
        Calculates the orientation of the GML surface based on its coordinates.
    
    unit_normal(a, b, c):
        Calculates the unit normal vector of a surface described by three points.
    """

    def __init__(self, gml_surface, surface_type, boundary=None):
        """
        Initializes the SurfaceGML object with GML coordinates and calculates its properties.

        Parameters
        ----------
        gml_surface : list
            A list of coordinates defining the surface. Should be a flat list of numbers in the order [x1, y1, z1, x2, y2, z2, ...].
        surface_type : str
            Type of surface element as str. (eg. WallSurface, RoofSurface, GroundSurface, etc...)
        boundary : str, optional
            The name or boundary identifier of the surface.
        """
        self.gml_surface = gml_surface
        self.surface_type = surface_type
        self.name = boundary
        self.surface_area = None
        self.surface_orientation = None
        self.surface_tilt = None

        try:
            self.surface_area = self.calculate_area()
            self.surface_orientation = self.get_gml_orientation()
            self.surface_tilt = self.get_gml_tilt()
        except Exception as e:
            print(f"Error processing surface '{self.name}': {e}")

    def fit_plane(self, points):
        """
        Fits a plane to the given 3D points and returns the normal vector and centroid.

        Parameters
        ----------
        points : np.ndarray
            An array of 3D points.

        Returns
        -------
        tuple
            A tuple containing the normal vector of the fitted plane and the centroid of the points.
        """
        centroid = points.mean(axis=0)
        centered_points = points - centroid
        u, s, vh = LA.svd(centered_points)
        normal = vh[-1, :]  # Use the last singular vector as the normal

        if LA.norm(normal) == 0:
            raise ValueError("Degenerate normal vector.")

        normal /= LA.norm(normal)

        return normal, centroid

    def project_to_2D(self, points, normal, centroid):
        """
        Projects the 3D points onto a 2D plane defined by the normal vector and centroid.

        Parameters
        ----------
        points : np.ndarray
            An array of 3D points.
        normal : np.ndarray
            The normal vector for the plane.
        centroid : np.ndarray
            The centroid of the given points.

        Returns
        -------
        np.ndarray or None
            The projected 2D coordinates or None if projection fails.
        """
        try:
            # Define the in-plane axes
            if not np.allclose(normal, [0, 0, 1]):
                in_plane_x = np.cross(normal, [0, 0, 1])
                if LA.norm(in_plane_x) == 0:
                    in_plane_x = np.array([1, 0, 0])  # Fallback if normal is parallel to z-axis
                else:
                    in_plane_x /= LA.norm(in_plane_x)
            else:
                in_plane_x = np.array([1, 0, 0])

            in_plane_y = np.cross(normal, in_plane_x)
            in_plane_y /= LA.norm(in_plane_y)

            # Translate points to the centroid
            translated = points - centroid

            # Project points onto the in-plane axes
            x_coords = np.dot(translated, in_plane_x)
            y_coords = np.dot(translated, in_plane_y)

            return np.vstack((x_coords, y_coords)).T
        except Exception as e:
            print(f"Error projecting points to plane: {e}")
            return None

    def calculate_area(self):
        """
        Calculates the area of the surface represented by the GML data.

        Returns
        -------
        float
            The calculated area of the surface.
        """
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

    def get_gml_tilt(self):
        """
        Calculates the tilt of the GML surface defined by its coordinates.

        The surface needs to be planar.

        Returns
        -------
        float
            The calculated tilt of the surface in degrees.
        """
        points = np.array(self.gml_surface).reshape(-1, 3)
        if len(points) < 3:
            print(f"Surface '{self.name}' has fewer than 3 points. Tilt set to 0.0.")
            return 0.0

        try:
            normal, _ = self.fit_plane(points)
            z_axis = np.array([0, 0, 1])

            tilt_rad = np.arccos(np.dot(normal, z_axis))
            tilt = np.degrees(tilt_rad)

            # Handle edge cases
            if np.isclose(tilt, 180.0):
                tilt = 0.0
            elif np.isnan(tilt):
                tilt = 90.0

            return tilt
        except Exception as e:
            print(f"Error calculating tilt for surface '{self.name}': {e}")
            return 0.0

    def get_gml_orientation(self):
        """
        Calculates the orientation of the GML surface based on its coordinates.

        The orientation returned is in TEASER coordinates.

        Returns
        -------
        float
            The calculated orientation of the surface.
        """
        points = np.array(self.gml_surface).reshape(-1, 3)
        num_points = len(points)

        if num_points < 3:
            print(f"Surface '{self.name}' has fewer than 3 points. Orientation set to None.")
            return None

        try:
            normal, _ = self.fit_plane(points)
            phi = np.degrees(np.arctan2(normal[1], normal[0]))
            #print(f"phi before: {phi}")
            if phi < 0:
                phi += 360
            teaser_orientation = phi
            #print(f"phi after: {phi}")
            #teaser_orientation = (90 - phi) % 360
            #print(f"orientation before: {teaser_orientation}")
            #if teaser_orientation > 180:
            #    teaser_orientation -= 360
            #print(f"orientation after: {teaser_orientation}")

            # ToDo: cleanup code. teaser conversion not necasary here. conversion is done in data/output/modelictemplate/conversion

            # Adjust orientation based on normal vector's z-component
            if np.isclose(abs(normal[2]), 1.0, atol=1e-3):
                
                if self.surface_type is not None:
                    if self.surface_type == 'RoofSurface':
                        orientation = -1  # Roof
                        return orientation
                    elif self.surface_type == 'GroundSurface':
                        orientation = -2  # Floor
                        return orientation
                    
                elif np.isclose(normal[2], 1.0, atol=1e-3):
                    orientation = -1  # Roof
                    return orientation
                
                elif np.isclose(normal[2], -1.0, atol=1e-3):
                    orientation = -2  # Floor
                    return orientation
            else:
                return teaser_orientation
        except Exception as e:
            print(f"Error calculating orientation for surface '{self.name}': {e}")
            return None

    def unit_normal(self, a, b, c):
        """
        Calculates the unit normal vector of a surface described by three points.

        Parameters
        ----------
        a : array-like
            Coordinates of the first point (x, y, z).
        b : array-like
            Coordinates of the second point (x, y, z).
        c : array-like
            Coordinates of the third point (x, y, z).

        Returns
        -------
        tuple
            A tuple representing the unit normal vector (x, y, z).
        """
        try:
            # Create vectors from the points
            ab = np.array(b) - np.array(a)
            ac = np.array(c) - np.array(a)

            # Calculate the cross product
            cross_prod = np.cross(ab, ac)

            # Calculate the magnitude of the cross product
            magnitude = LA.norm(cross_prod)
            if magnitude == 0:
                print("Degenerate normal vector. Returning (0, 0, 0).")
                return (0.0, 0.0, 0.0)

            # Return the unit normal vector
            return tuple(cross_prod / magnitude)
        except Exception as e:
            print(f"Error calculating unit normal: {e}")
            return (0.0, 0.0, 0.0)