
import numpy as np
from mmesh import Mesh
from slfpy import SLF
from shapely.geometry import Point,Polygon
from mshapely import DF

import matplotlib.pyplot as plt

def test_io():
  slf=SLF()
  slf.writeSLF("../data/slf/example1.slf")
  _m=Mesh("../data/slf/example1.slf")
  fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
  fig.tight_layout()
  _m.boundaries.plot(axe=axes).savePlot("../data/mesh/boundaries.1.png")

def test_boundaries():
  
  polygon = Point((0,0)).buffer(100)
  hole1 = Point((-50,0)).buffer(20)
  hole2 = Point((50,0)).buffer(20)
  polygon = Polygon(polygon.exterior,[hole1.exterior.coords[::-1],hole2.exterior.coords[::-1]])
  distance=polygon.inearest(maxDistance=100,angle=90)
  distance[:,2]=DF.getD_l(1,1.2,distance[:,2])
  density = np.column_stack((distance,np.ones(len(distance))*1.2))
  
  df=DF(density,minDensity=1,maxDensity=100)
  
  polygon=polygon.dresample(df)
  
  distance=polygon.inearest(maxDistance=100,angle=90)
  distance[:,2]=DF.getD_l(1,1.2,distance[:,2],1)
  
  density = np.column_stack((distance,np.ones(len(distance))*1.2))
  df=DF(density,minDensity=1,maxDensity=100)
  polygon.msh("test/data/test.msh",df)
  
  fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
  fig.tight_layout()
  _m=Mesh("test/data/test.msh")
  _m.plotBoundaries(axe=axes).plot(axe=axes).savePlot("test/data/test.mesh.2.png")

def test_edges_geojson():
  _m=Mesh("test/data/test.msh")
  _m.writeEdges("test/data/edge.geojson")
  _m.writeBoundaries("test/data/boundary.geojson")
  
    

if __name__ == "__main__":
  # test_io()
  # test_boundaries()
  test_edges_geojson()