import os
import numpy as np
from matplotlib.tri import Triangulation
from shapely.geometry import Point,LineString,Polygon,MultiPoint,MultiLineString,MultiPolygon,GeometryCollection
import mshapely
from mshapely.misc import add_method

from .io import createGEO,createMSH

#
# Create Gmsh
#
@add_method(GeometryCollection)
def msh(self,*args,**kwargs):
  return self.toShape().msh(*args,**kwargs)


@add_method([Point,MultiPoint,LineString,MultiLineString])
def msh(self,*args,**kwargs):
  return self
  
@add_method(Polygon)
def msh(self,path,density,*args,**kwargs):
  geo = createGEO(self,path,density,*args,**kwargs)
  return createMSH(geo,path)

@add_method(MultiPolygon)
def msh(self,*args,**kwargs):
  geo=self.largest()
  geo.msh(*args,**kwargs)
  return geo
  

