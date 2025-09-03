import pyvista as pyv
import numpy as np
import math
import vtk
import os
from dataclasses import dataclass
from typing import Union, Tuple

# Generating spiral through a sphere https://perswww.kuleuven.be/~u0017946/publications/Papers97/art97a-Saff-Kuijlaars-MI/Saff-Kuijlaars-MathIntel97.pdf
def create_sphere_points(center, radius, num_points):
    # Generate points along a spiral trajectory on the surface of the sphere
    inv_sqrt_n = 1.0 / math.sqrt(num_points)
    phi = 0.0
    points = []

    for i in range(num_points):
        # Compute polar coordinates of marker positions
        ck = -1.0 + (2.0 * i) / (num_points - 1)
        theta = math.acos(ck)

        if i == 0 or i == num_points - 1:
            phi = 0.0
        else:
            phi = (phi + 3.6 * inv_sqrt_n / math.sqrt(1 - ck**2)) % (2 * math.pi)

        # Spiral around z. Convert to Cartesian coordinates
        # x = center[0] + radius * np.sin(theta) * np.cos(phi)
        # y = center[1] + radius * np.sin(theta) * np.sin(phi)
        # z = center[2] + radius * np.cos(theta)

        # Spiral around x. Convert to Cartesian coordinates
        y = center[1] + radius * np.sin(theta) * np.cos(phi)
        z = center[2] + radius * np.sin(theta) * np.sin(phi)
        x = center[0] + radius * np.cos(theta)

        # Append the coordinates to the points list
        points.append([x, y, z,])

    # Convert the list to a NumPy array
    points = np.array(points)
    return points



# Create two point clouds representing spheres with different radii
center = (36, 40, 10)
radius2 = 4.0 # 5um in cm outer radius
radius1 = 4.0 - 10.0/20  # make sure the thickness equal to marker distance, 0.5um, on flagellum
num_points = 48

points1 = create_sphere_points(center, radius1, num_points)
points2 = create_sphere_points(center, radius2, num_points)
point=np.vstack((points1,points2))
#point_all=np.vstack((points1,points2,marker_points))
#print("coordinates:")
#print("1:",points1)
#print("2:",points2)
#print("3:",point)
print("Fenix was here :)")

# Create two PyVista point clouds

point_cloud = pyv.PolyData(point)

# Create arrays for IDs and add them as scalars to the point clouds for color-coding

ids1 = np.arange(len(points1))
ids2 = np.arange(len(points1), len(point))
print(ids1, ids2)
point_cloud['Point ID'] = np.concatenate((ids1, ids2))



# Compute the 3D Delaunay triangulation for each point cloud
delaunay = point_cloud.delaunay_3d(alpha=36.0,tol=0.000001,offset=2.5,progress_bar=False)
edges=delaunay.extract_all_edges()

#get neighbouring particles,id and coordinates
#This creates a dictionary that can be called.
@dataclass
class Point:
    id: int
    coord: Union[None, Tuple[float]]
    neighbors: Tuple[int]

def map_points(delaunay):
    pt_map: Dict[int, Point] = dict()
    for i in range(delaunay.n_points):
        pt_map[i] = Point(
            id=i,
            coord=tuple(delaunay.points[i, :]),
            neighbors=tuple(delaunay.point_neighbors(i))

        )
    return pt_map
pts = map_points(delaunay)


for i in range(delaunay.n_cells):
    cell = delaunay.get_cell(i)
    for j in range(cell.GetNumberOfPoints()):
        coord = cell.GetPoints().GetPoint(j)
        idx = cell.point_ids[j]


#Generates particles.dat

scaler=1.0e-04  #for scaling any coord,bond length etc from um to cm
move=0; # ??? 1e-03 # move the cell body

with open("particles.dat", "w") as f:
    for idx in pts.keys():
        i = pts[idx].id
        p0 = np.array(pts[idx].coord)
        f.write(f"{i}\t{(p0[0]*scaler)+move}\t{(p0[1]*scaler+move)}\t{(p0[2]*scaler)+move}\t1\t0\t1\n")


if os.path.exists("particles.dat"):
    print("Particles.dat created!!!!")
else:
    print("Failed to create particles.dat!!!! :(")




#For finding the max number of neighbour of point have.
max_nn = 0
for _, pt in pts.items():
    max_nn = max(max_nn, len(pt.neighbors))



#Generates bonds.csv

with open("bonds.csv","w") as f:
    k=4e-4/max_nn
    for idx in pts.keys():
        i = pts[idx].id
        p0 = np.array(pts[idx].coord)
        for j in pts[idx].neighbors:
            p1 = np.array(pts[j].coord)
            dist = np.sqrt(np.sum(((p1-p0)*scaler)**2))
            f.write(f"{i}\t{j}\t{k}\t{dist}\n")
if os.path.exists("particles.dat"):
    print("bonds.csv created!!!!")
else:
    print("Failed to create bonds.csv!!! :(")


marker_points = []
y = 10*scaler
z = 10*scaler
for i in range(23):
    # Compute x of marker positions
    x = (20 + 0.5 * i)*scaler
    marker_points.append([x, y, z,])
marker_points = np.array(marker_points)


# Plot the point clouds and their Delaunay triangulations using PyVista's plotting capabilities
edges.plot(line_width=1,color='k')
p = pyv.Plotter()
p.add_axes()
p.add_mesh(delaunay, color='lightblue', opacity=0.6, show_edges=True)
p.add_points(point_cloud,scalars='Point ID', cmap='coolwarm', point_size=6)
p.add_mesh(pyv.Sphere(center=center, radius=radius1), color='red', opacity=0.2)
p.add_mesh(pyv.Sphere(center=center, radius=radius2), color='green', opacity=0.2)
p.show()


p1 = pyv.Plotter()
#p1.add_points(point_cloud, color = 'red', point_size=6)
p1.add_points(marker_points, color = 'green', point_size=6)
p1.add_axes()
p1.show()


#extra stuff for optimization or settings
#print(f"{i}\t{j}\t{k/len(pts[idx].neighbors)}\t{dist}")   # this equation will print the neighbouring particles and even it by it's number to provide optimal k etc. id_point=1 neighbour_points=2 so k/n_point

#For finding the max number of neighbour of point have.
#max_nn = 0
#for _, pt in pts.items():
#    max_nn = max(max_nn, len(pt.neighbors))







