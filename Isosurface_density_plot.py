import numpy as np
from scipy import stats
from mayavi import mlab
import pandas as pd


#Select csv files as input (dfn) the data column to be plotted (x,y,z)
df1 = pd.read_csv()
df2 = pd.read_csv()
x = df1['Center of Mass (Bounding box): X (px)']
x = np.append(x, df2['Center of Mass (Bounding box): X (px)'])
xmax = np.amax(x)
xmin = np.amin(x)
x = (x-xmin)/(xmax-xmin)
y = df1['Center of Mass (Bounding box): Y (px)']
y = np.append(y, df2['Center of Mass (Bounding box): Y (px)'])
ymax = np.amax(y)
ymin = np.amin(y)
y = (y-ymin)/(ymax-ymin)
z = df1['Center of Mass (Bounding box): Z (px)']
z = np.append(z, df2['Center of Mass (Bounding box): Z (px)'])
zmax = np.amax(z)
zmin = np.amin(z)
z = (z-zmin)/(zmax-zmin)

values = np.array([x, y, z])

kde = stats.gaussian_kde(values)

xmin, ymin, zmin = values.T.min(axis=0)
xmax, ymax, zmax = values.T.max(axis=0)

xi, yi, zi = np.mgrid[xmin:xmax:70j, ymin:ymax:70j, zmin:zmax:70j]

coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
density = kde(coords).reshape(xi.shape)

print (np.amin(density), np.amax(density))


#maxid = np.where(density == np.amax(density))
#maxid = np.argmax(density)
#print (coords.shape, maxid)

#Characterstics of contour plot
mlab.contour3d(xi, yi, zi, density, contours=20, vmin=0, vmax=15,
               opacity=0.4)
mlab.axes()
mlab.show()

