import os
import numpy as np
from mmesh import MESH,MSLF,MMSH
from slfpy import SLF
from shapely.geometry import Point,Polygon
from mshapely import DF

import matplotlib.pyplot as plt

def test_MESH():
  grid=SLF.createGrid()
  mesh=MESH(grid['xy'],grid['ikle'])
  mesh.write("slf.1.slf")
  mesh.write("mesh.1.msh")
  mesh1=MSLF("slf.1.slf")
  mesh2=MMSH("mesh.1.msh")
  
  np.testing.assert_almost_equal(mesh.x,mesh1.x)
  np.testing.assert_almost_equal(mesh.y,mesh1.y)
  np.testing.assert_almost_equal(mesh.triangles,mesh1.triangles)
  
  np.testing.assert_almost_equal(mesh.x,mesh2.x)
  np.testing.assert_almost_equal(mesh.y,mesh2.y)
  np.testing.assert_almost_equal(mesh.triangles,mesh2.triangles)
  os.remove("slf.1.slf")
  os.remove("mesh.1.msh")
  
if __name__ == "__main__":
  test_MESH()