import os
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from shapely.geometry import Point,LineString,Polygon,MultiPoint,MultiLineString,MultiPolygon,GeometryCollection
import mshapely
from mshapely.plot import plotTri,plotSave

from slfpy import SLF
from .gmsh import GMSH,Mesh

class MESH(Triangulation):
  """
  Mesh object for MERCAN.
  It's an extension of matplotlib.tri.Triangulation.
  It contains new attributes and methods such as writing boundaries and edges to shapely.
  
  Note
  ----
  x:x-coordinate
  y:y-coordinate
  n1:point/node id-1 of triangle/element
  n2:point/node id-2 of triangle/element
  n3:point/node id-3 of triangle/element
  
  Parameters
  ----------
  xy:2D ndarray, [[x,y]]
  triangles:2D ndarray, [[n1,n2,n3]]
  
  """
  def __init__(self,xy,*args,**kwargs):
    x=xy[:,0]
    y=xy[:,1]
    super().__init__(x,y,*args,**kwargs)
    self._boundaries=None
  
  @property
  def geoedges(self):
    """
    Converts numpy edges to LineString object
    """
    xy=np.column_stack((self.x,self.y))
    edges=xy[self.edges]
    return GeometryCollection(list(map(LineString,edges)))
    
  
  @property
  def boundaries(self):
    """
    Gets boundaries/outline of the mesh and converts them in a Polygon object
    """
    if self._boundaries is None:
      neighbors = self.get_cpp_triangulation().get_neighbors()
      
      a=np.column_stack((self.triangles[neighbors[:,0]<0,0],self.triangles[neighbors[:,0]<0,1]))
      b=np.column_stack((self.triangles[neighbors[:,1]<0,1],self.triangles[neighbors[:,1]<0,2]))
      c=np.column_stack((self.triangles[neighbors[:,2]<0,2],self.triangles[neighbors[:,2]<0,0]))
      bedges=np.concatenate((a,b,c))
      
      pols=[]
      pol=[]
      index=-1;start=-1;end=-1
      while(len(bedges)>0):
        index=np.where(bedges==end)[0];
        if(len(index)==0):
          if(len(pol)>0):
            pols.append(pol)
            pol=[]
          start=bedges[0][0];
          end=bedges[0][1];
          pol.append([self.x[start],self.y[start]])
          pol.append([self.x[end],self.y[end]])
          bedges=np.delete(bedges, 0, 0)
        else:
          index=index[0]
          end= bedges[index][1] if bedges[index][0]==end else bedges[index][0]
          pol.append([self.x[end],self.y[end]]);
          bedges=np.delete(bedges, index, 0)
          if len(bedges)==0 and len(pol)>0:pols.append(pol)
        
      multi=MultiPolygon([Polygon(pol) for pol in pols])
      boundaries,holes = multi.largest(return_other=True)
      if len(holes)!=0:boundaries=boundaries.difference(MultiPolygon(holes))
      self._boundaries = boundaries
    return self._boundaries
    
  @property  
  def xy(self):
    """
    Get xy array
    """
    return np.column_stack((self.x,self.y))
  
  def setXY(self,x,y):
    """
    Set x and y coordinates
    Parameters
    ----------
    x:1D ndarray, [x]
    y:1D ndarray, [y]
    """
    self.x=x
    self.y=y
    return self
      
  def write(self,path,**kwargs):
    """
    Write MESH to file depedning on the file extension.
    SELAFIN and MSH are acceptable
    
    Parameters
    ----------
    path:str
      Extention must have (.slf) or (.msh)
    
    """
    ext = os.path.splitext(path)[1]
    if ext==".slf":
      SLF().addMesh(self.xy,self.triangles,**kwargs).write(path)
      return self
    elif ext==".msh":
      mesh=Mesh(self.xy,[("triangle", self.triangles)])
      GMSH.write(path,mesh)
      return self
      
    raise Exception("Method does not exist for {}".format(ext))
  
  @staticmethod
  def delete(path):
    """
    Delete file
    """
    if os.path.exists(path):os.remove(path)
    
  def plot(self,**kwargs):
    """Plot MESH with matplotlib
    """
    plotTri(self,**kwargs)
    return self  

  def plotBoundaries(self,*args,**kwargs):
    """Plot boundaries with matplotlib
    """    
    self.boundaries.plot(*args,**kwargs)
    return self

  def savePlot(self,path):
    """Save plot to file
    """
    plotSave(path)
    return self

  @staticmethod
  def read(path):
    """
    Read MESH file using file extention
    """
    ext = os.path.splitext(path)[1]
    if ext==".slf":return MSLF(path)
    if ext==".msh":return MMSH(path)
    raise Exception("Format does not exist")

class MSLF(MESH,SLF):
  """
  Hybrid object of MESH object and SLF object
  """
  def __init__(self,path):
    SLF.__init__(self,path)
    MESH.__init__(self, np.column_stack((self.MESHX,self.MESHY)),self.IKLE2)
    
  # TODO: write values  
  # def write(self,path,**kwargs):
  #   ext = os.path.splitext(path)[1]
  #   if ext==".msh":return MESH.write(self,path)
  #   elif ext==".slf":return SLF.write(self,path)
      
  #   raise Exception("Method does not exist for {}".format(ext))
  
class MMSH(MESH):
  """
  Hybrid object of MESH object and MSH object
  """  
  def __init__(self,path):
    msh=GMSH.read(path)
    MESH.__init__(self, msh.points,msh.get_cells_type("triangle"))