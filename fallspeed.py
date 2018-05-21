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
from glob import glob
#import gzip

#shapefile='subsequent-0.1/dendrite-650e-6-720-0.1-subsequent-7cVE9G0d.agg'
#metadata='subsequent-0.1/dendrite-650e-6-720-0.1-subsequent-7cVE9G0d.agg.gz.meta'

#shapefile='./simultaneous-0.0/dendrite-650e-6-18-0.0-simultaneous-p5stn3Dr.agg'
#metadata='./simultaneous-0.0/dendrite-650e-6-18-0.0-simultaneous-p5stn3Dr.agg.gz.meta'

#shapefile='./rimeonly_120e-6/dendrite-650e-6-1-0.1-rimeonly_120e-6-YWrwsxzK.agg'
#metadata='./rimeonly_120e-6/dendrite-650e-6-1-0.1-rimeonly_120e-6-YWrwsxzK.agg.gz.meta'

xy=(0,1)
yz=(1,2)
xz=(0,2)
cols=['shapefile','Jdmax','Ddmax','Dprojcirc','projarea','mass','vDJ','vDD','vDproj']

path = '/data/optimice/scattering_databases/shape_files_jussi/' # './'
folders=['simultaneous-0.0/','simultaneous-0.1/','simultaneous-0.2/','simultaneous-0.5/','simultaneous-1.0/','simultaneous-2.0/','subsequent-0.1/','subsequent-0.2/','subsequent-0.5/','subsequent-1.0/','subsequent-2.0/','rimeonly_120e-6']

for folder in folders:
    shpfiles = glob(path+folder+'*.agg')
    data = pd.DataFrame(columns=cols)
    jj = 0
    for shapefile in shpfiles:
        shape=np.loadtxt(shapefile)
        metadata=shapefile+'.gz.meta'
        attributes=eval(open(metadata).read())
        d=attributes['grid_res']
        try:
            hull3d=ConvexHull(shape)
            hull3d=hull3d.points[hull3d.vertices]
        except:
            hull3d=shape
        dmax = 0
        for i in range(0,hull3d.shape[0]-1):
            p0 = hull3d[i,:]
            for j in range(i+1,hull3d.shape[0]):
                p1 = hull3d[j,:]
                r=p0-p1
                dist = np.dot(r,r)
                if dist > dmax:
                    dmax = dist
        dmax = d*dmax**0.5

        try:
            projection=pd.DataFrame(shape[:,xy]).drop_duplicates().values
            hull=ConvexHull(projection)
            circle=make_circle(hull.points[hull.vertices])
        except:
            circle=make_circle(projection)

#    plt.figure()
#    ax=plt.gca()
#    ax.scatter(projection[:,0],projection[:,1])
#    pltcircle=plt.Circle(circle[0:2],circle[2],alpha=0.2,color='g')
#    ax.add_artist(pltcircle)
#    plt.axis('equal')

        area=projection.shape[0]*d**2.0
        diam=d*circle[2]*2.0
        dmax_att=attributes['max_diam']
        area_ratio=projection.shape[0]*d**2.0/(np.pi*(d*circle[2])**2.0)
        mass=shape.shape[0]*d**3.0*916.0

        from heymsfield_2010 import dia2vel
        v=dia2vel(diam, 1.0, 1.6e-5, mass, area)
        print(folder[-4:],jj,len(shpfiles),v,dia2vel(dmax, 1.0, 1.6e-5, mass, area),dia2vel(attributes['max_diam'], 1.0, 1.6e-5, mass, area))
        data.loc[jj]=[shapefile,dmax_att,dmax,diam,area,mass,v,dia2vel(dmax, 1.0, 1.6e-5, mass, area),dia2vel(attributes['max_diam'], 1.0, 1.6e-5, mass, area)]
        jj = jj + 1
    data.to_csv(folder[2:-1]+'.csv')
