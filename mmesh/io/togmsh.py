import numpy as np
import subprocess
import os
# from ..gmsh import GMSH
from ..mesh import MMSH
from mshapely import DF

def createMSH(input,output):
  # command = "gmsh {0} -2 -format msh2 -algo frontal -smooth {2} -o {1}.msh".format(file,name,smooth)
  input = os.path.splitext(input)[0]+".geo"
  output = os.path.splitext(output)[0]+".msh"
  command = "gmsh {0} -2 -smooth 10 -algo frontal -o  {1}".format(input,output)
  print(command)
  subprocess.call(command, shell=True)
  return MMSH(output)
    
def createGEO(geo,path,df,*args,**kwargs):
  output = os.path.splitext(path)[0]+".geo"
  points=geo.np
  strPoints = ""
  strLines = ""
  strLoop = ""
  strSurface = ""
  PID = 1
  LID = 1
  LLID = 1
  
  # Create density points that are not in shoreline
  if(np.any(np.unique(df.dp[:,4])>1)):
    dp = df.dp
    xy=dp[dp[:,4]==0,:2]
    for p0 in xy:
      strPoints += "Point({0:n}) = {{{1},{2},{3},{4}}};\n".format(PID, p0[0], p0[1], 0, 100)
      PID +=1
    
  for ipolygon in np.unique(points[:, 0]):
    subpoints = points[points[:, 0] == ipolygon]
    
    istart=LID
    IPID=PID
    for i in range(len(subpoints)-1):
      p0=subpoints[i]
      _TID= PID+1 if i < len(subpoints)-2 else IPID
      # if(PID==32):print(p0)
      strPoints += "Point({0:n}) = {{{1},{2},{3},{4}}};\n".format(PID, p0[2], p0[3], 0, 100)
      strLines += "Line({0:n}) = {{{1:n},{2:n}}};\n".format(LID, PID, _TID)
      PID +=1
      LID +=1

    iend = LID - 1
    strLoop += "l{0:n} = newreg; ".format(LLID)
    strLoop += "Line Loop(l{0:n}) = {{{1:n}:{2:n}}};\n".format(LLID, istart, iend)
    LLID +=1
  
  strTs = ""
  for n in range(1, LLID):
    strTs += "l{0:n},".format(n)
  strTs = strTs[:-1]
  strSurface +="s{0:n} = newreg;".format(1)
  strSurface +="Plane Surface(s{0}) = {{{1}}};".format(1, strTs)
  
  with open(output,"w") as geofile:
    geofile.write("{0}\n".format(strPoints))
    geofile.write("{0}\n".format(strLines))
    geofile.write("{0}\n".format(strLoop))
    geofile.write("{0}\n".format(strSurface))
    geofile.write("{0}\n".format(getAttractors(points,df)))
  
  
  return output

  
def getAttractors(points,df):
  
  
  minDensity=df.minDensity
  maxDensity=df.maxDensity
  strAttractor =""
  ATID=1
  gATID=[str(1)]
  dp=df.dp
  density=dp[:,2]
  growth=dp[:,3]
  distances=DF.getl_D(density,growth,maxDensity)
  
  # Get # of density points that are not in shoreline
  nxy=len(dp[dp[:,4]==0])
  
  for i,point in enumerate(df.dp):
    ATID +=1 
    pointID = point[5] if point[4]==0 else nxy+point[5] # Check if density is in shoreline or not
    
    strAttractor +="Field[{0}] = Attractor;Field[{0}].NodesList = {{{1:n}}};".format(ATID,pointID+1)    

    distance=distances[i]
    _d=density[i]
    ATID +=1
    strAttractor +="Field[{0}] = Threshold;Field[{0}].IField = {1};Field[{0}].DistMax = {4};Field[{0}].DistMin = 0;Field[{0}].LcMax = {3};Field[{0}].LcMin = {2:.1f};\n".format(ATID,ATID-1,_d,maxDensity,distance)
    gATID.append(str(ATID))
  
  ATID += 1
  strAttractor +='Field[1] = MathEval;Field[1].F = "{}";'.format(maxDensity)
  strAttractor +="Field[{0}] = Min;Field[{0}].FieldsList = {{{1}}};\n".format(ATID,",".join(gATID))
  strAttractor +="Background Field = {0};\n".format(ATID)
  strAttractor +="Mesh.LcIntegrationPrecision = 1e-3;\n"
  strAttractor +="Mesh.CharacteristicLengthExtendFromBoundary = 0;\n"
  strAttractor +="Mesh.CharacteristicLengthFromPoints = 0;\n"
  
  return strAttractor
  
  
# def getAttractors(points):
#   attractors=np.concatenate([
#     np.arange(1E1,1E2,1E1),
#     np.arange(1E2,1E3,1E2),
#     np.arange(1E3, 1E4, 1E3),
#     np.arange(1E4, 1E5, 1E4),
#     ])
    
#   # pp=attractors[np.abs(attractors[:,None]-points).argmin(axis=0)]
#   index=np.abs(attractors[:,None]-points).argmin(axis=0)
  
#   strAttractorTreshold =""
#   ATID=0
#   for i,density in enumerate(attractors):
#     nodelist=np.where(i==index)[0]
#     if len(nodelist)>0:
#       strnodelist = ",".join(nodelist.astype('str'))
#       ATID +=1 
#       strAttractorTreshold +="Field[{0}] = Attractor;Field[{0}].NodesList = {{{1}}};".format(ATID,strnodelist)
#       n= np.maximum(np.floor(np.log(maxDensity/density)/np.log(growth)-1),1)
#       maxDistance = (minDensity*np.power(growth,n+1)-minDensity)/(growth-1)  
#       ATID +=1
#       strAttractorTreshold +="Field[{0}] = Threshold;Field[{0}].IField = {1};Field[{0}].DistMax = {4};Field[{0}].DistMin = 0;Field[{0}].LcMax = {3};Field[{0}].LcMin = {2:.1f};\n".format(ATID,ATID-1,density,maxDensity,distance)
#       print(strnodelist)
    
    
#     # print(index[i==index])
#   # print(np.abs(attractors[:,None]-points).argmin(axis=0))
  
  None