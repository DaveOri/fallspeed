# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:02:37 2018

@author: dori
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from small_circle import make_circle

shapefile='subsequent-0.1/dendrite-650e-6-720-0.1-subsequent-7cVE9G0d.agg'
metadata='subsequent-0.1/dendrite-650e-6-720-0.1-subsequent-7cVE9G0d.agg.gz.meta'

shape=np.loadtxt(shapefile)
attributes=eval(open(metadata).read())
d=attributes['grid_res']
hull3d=ConvexHull(shape)
hull3d=hull3d.points[hull3d.vertices]

dmax = 0
for i in range(1,hull3d.shape[0]):
    p1 = hull3d[i,:]
    p0 = hull3d[i-1,:]
    dist = np.dot(p0,p1)
    if dist > d:
        dmax = dist
dmax = d*dmax**0.5

xy=(0,1)
yz=(1,2)
xz=(0,2)
projection=pd.DataFrame(shape[:,xz]).drop_duplicates().values
hull=ConvexHull(projection)
circle=make_circle(hull.points[hull.vertices])

plt.figure()
ax=plt.gca()
ax.scatter(projection[:,0],projection[:,1])
pltcircle=plt.Circle(circle[0:2],circle[2],alpha=0.2,color='g')
ax.add_artist(pltcircle)
plt.axis('equal')

area=projection.shape[0]*d**2.0
diam=d*circle[2]*0.5
#diam=attributes['max_diam']
area_ratio=projection.shape[0]*d**2.0/(np.pi*(d*circle[2])**2.0)
mass=shape.shape[0]*d**3.0*916.0
from heymsfield_2010 import dia2vel

v=dia2vel(diam, 1.0, 1.6e-5, mass, area)
print(v)
print(dia2vel(dmax, 1.0, 1.6e-5, mass, area))
print(dia2vel(attributes['max_diam'], 1.0, 1.6e-5, mass, area))