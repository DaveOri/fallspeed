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

xy=(0,1)
yz=(1,2)
xz=(0,2)
projection=pd.DataFrame(shape[:,xy]).drop_duplicates().values
hull=ConvexHull(projection)
circle=make_circle(hull.points[hull.vertices])

plt.figure()
ax=plt.gca()
ax.scatter(projection[:,0],projection[:,1])
pltcircle=plt.Circle(circle[0:2],circle[2],alpha=0.2,color='g')
ax.add_artist(pltcircle)
plt.axis('equal')

d=attributes['grid_res']
area_ratio=projection.shape[0]*d**2.0/(np.pi*(d*circle[2])**2.0)