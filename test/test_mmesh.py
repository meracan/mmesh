
import numpy as np
from mmesh import MESH,MSLF,MMSH
from slfpy import SLF
from shapely.geometry import Point,Polygon
from mshapely import DF

import matplotlib.pyplot as plt

def test_mesh():
  grid=SLF.createGrid()
  mesh=MESH(grid['xy'],grid['ikle'])
  mesh.boundaries.write("test/data/boundaries.1.geojson")
  mesh.geoedges.write("test/data/edges.1.geojson")
  mesh.geoedges.plot().savePlot("test/data/mesh.1.png")
  mesh.boundaries.plot(showPoints=True,pointStyle={"c":"black"}).savePlot("test/data/boundaries.1.png")
  mesh.write("test/data/slf.1.slf")
  mesh.write("test/data/mesh.1.msh")
  
  mesh1=MSLF("test/data/slf.1.slf")
  mesh2=MMSH("test/data/mesh.1.msh")
  
  np.testing.assert_almost_equal(mesh.x,mesh1.x)
  np.testing.assert_almost_equal(mesh.y,mesh1.y)
  np.testing.assert_almost_equal(mesh.triangles,mesh1.triangles)
  
  np.testing.assert_almost_equal(mesh.x,mesh2.x)
  np.testing.assert_almost_equal(mesh.y,mesh2.y)
  np.testing.assert_almost_equal(mesh.triangles,mesh2.triangles)
  
  
  


if __name__ == "__main__":
  test_mesh()
  # test_mslf()
  # test_mmesh()