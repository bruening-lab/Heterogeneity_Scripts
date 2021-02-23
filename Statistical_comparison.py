
import numpy as np
from scipy import stats
from mayavi import mlab
import pandas as pd
import scipy.ndimage as sci
import matplotlib.pyplot as plt

#Select csv files as input 
GroupA = []

GroupB = []

xmax_ = -1000
xmin_ = 1000
ymax_ = -1000
ymin_ = 1000
zmax_ = -1000
zmin_ = 1000

#For GroupA select the column of data to plot (x,y,z)
all1 = [[], [], [], [], [], [], [], []]
allex1 = [[], [], [], [], [], [], [], []]
cc = 0
l1 = 0.0
for i in glp1r:
	df1 = pd.read_csv(i)
	x = df1['Center of Mass (Bounding box): X (px)']
	xmax = np.amax(x)
	if xmax > xmax_:
		xmax_ = xmax
	xmin = np.amin(x)
	if xmin < xmin_:
		xmin_ = xmin
	y = df1['Center of Mass (Bounding box): Y (px)']
	ymax = np.amax(y)
	if ymax > ymax_:
		ymax_ = ymax
	ymin = np.amin(y)
	if ymin < ymin_:
		ymin_ = ymin
	z = df1['Center of Mass (Bounding box): Z (px)']
	zmax = np.amax(z)
	if zmax > zmax_:
		zmax_ = zmax
	zmin = np.amin(z)
	if zmin < zmin_:
		zmin_ = zmin
	values = np.array([x, y, z])
	l1 += len(x)/np.float(len(glp1r))
	extrema = np.array([[xmin, ymin, zmin], [xmax, ymax, zmax]])
	all1[cc] = values
	allex1[cc] = extrema
	cc+=1

#For GroupB select the column of data to plot (x,y,z)
all2 = [[], [], [], [], [], [], []]
allex2 = [[], [], [], [], [], [], []]
cc = 0
l2 = 0.0
for i in lepr:
	df1 = pd.read_csv(i)
	x = df1['Center of Mass (Bounding box): X (px)']
	xmax = np.amax(x)
	if xmax > xmax_:
		xmax_ = xmax
	xmin = np.amin(x)
	if xmin < xmin_:
		xmin_ = xmin
	y = df1['Center of Mass (Bounding box): Y (px)']
	ymax = np.amax(y)
	if ymax > ymax_:
		ymax_ = ymax
	ymin = np.amin(y)
	if ymin < ymin_:
		ymin_ = ymin
	z = df1['Center of Mass (Bounding box): Z (px)']
	zmax = np.amax(z)
	if zmax > zmax_:
		zmax_ = zmax
	zmin = np.amin(z)
	if zmin < zmin_:
		zmin_ = zmin
	values = np.array([x, y, z])
	l2 += len(x)/np.float(len(lepr))
	extrema = np.array([[xmin, ymin, zmin], [xmax, ymax, zmax]])
	all2[cc] = values
	allex2[cc] = extrema
	cc+=1	

cut = 0.1
step = 20
stepp = 10
deltax = (xmax_ - xmin_)/step
deltay = (ymax_ - ymin_)/stepp
deltaz = (zmax_ - zmin_)/step
density1 = np.zeros((len(all1), step, stepp, step))
density2 = np.zeros((len(all2), step, stepp, step))
test = np.zeros((step, stepp, step))

polariy = -1 # +1 (glp1r > lepr), -1 (glp1r < lepr)

scale = l1 / l2
print 'scaling = mean(glp1r)/mean(lepr)'
print scale

pointx = []
pointy = []
pointz = []
pointt = []
for i in range(step):
	for j in range(stepp):
		for k in range(step):
			x = xmin_ + deltax*i
			y = ymin_ + deltay*j
			z = zmin_ + deltaz*k
			xx = xmin_ + deltax*(i+1)
			yy = ymin_ + deltay*(j+1)
			zz = zmin_ + deltaz*(k+1)	
			for c in range(len(all1)):
				#print (c)
				for cc in range(len(all1[c][0]) - 1):
					if (x <= all1[c][0, cc] and y <= all1[c][1, cc] and z <= all1[c][2, cc]):
						if (xx > all1[c][0, cc] and yy > all1[c][1, cc] and zz > all1[c][2, cc]):
							density1[c, i, j, k] += 1
			for c in range(len(all2)):
				#print (c)
				for cc in range(len(all2[c][0]) - 1):
					if (x <= all2[c][0, cc] and y <= all2[c][1, cc] and z <= all2[c][2, cc]):
						if (xx > all2[c][0, cc] and yy > all2[c][1, cc] and zz > all2[c][2, cc]):
							density2[c, i, j, k] += 1						
			dummy, t = stats.ttest_ind(1.0*density1[:, i, j, k], scale*density2[:, i, j, k], equal_var = False)
			if np.isnan(t) == False and t > 0.0:
				test[i, j, k] = t
				if t <= cut and polariy*(np.sum(1.0*density1[:, i, j, k]) - np.sum(scale*density2[:, i, j, k])) > 0.0:
					pointx.append((x+xx)/2.0)
					pointy.append((y+yy)/2.0)
					pointz.append((z+zz)/2.0)
					pointt.append(1.0-t)
			elif (np.sum(density1[:, i, j, k])+np.sum(density2[:, i, j, k]) > 0.0):
				test[i, j, k] = 1.0

#print pointx, pointy, pointz
test = sci.filters.gaussian_filter(test, (1.5, 1.5, 1.5));

mlab.points3d(pointx, pointy, pointz, pointt, mode='sphere', colormap='autumn') # mode = 'sphere'

polariy = +1 # +1 (glp1r > lepr), -1 (glp1r < lepr)

pointx = []
pointy = []
pointz = []
pointt = []
for i in range(step):
	for j in range(stepp):
		for k in range(step):
			x = xmin_ + deltax*i
			y = ymin_ + deltay*j
			z = zmin_ + deltaz*k
			xx = xmin_ + deltax*(i+1)
			yy = ymin_ + deltay*(j+1)
			zz = zmin_ + deltaz*(k+1)							
			dummy, t = stats.ttest_ind(1.0*density1[:, i, j, k], scale*density2[:, i, j, k], equal_var = False)
			if np.isnan(t) == False and t > 0.0:
				test[i, j, k] = t
				if t <= cut and polariy*(np.sum(1.0*density1[:, i, j, k]) - np.sum(scale*density2[:, i, j, k])) > 0.0:
					pointx.append((x+xx)/2.0)
					pointy.append((y+yy)/2.0)
					pointz.append((z+zz)/2.0)
					pointt.append(1.0-t)
			elif (np.sum(density1[:, i, j, k])+np.sum(density2[:, i, j, k]) > 0.0):
				test[i, j, k] = 1.0

#print pointx, pointy, pointz
test = sci.filters.gaussian_filter(test, (1.5, 1.5, 1.5));

xi, yi, zi = np.mgrid[xmin_:xmax_:20j, ymin_:ymax_:10j, zmin_:zmax_:20j]

mlab.contour3d(xi, yi, zi, test, contours=10, vmin=0.0, vmax=1.0, opacity=0.1)
mlab.points3d(pointx, pointy, pointz, pointt, mode='sphere', colormap='winter') # mode = 'cube'
mlab.axes()
mlab.view(40, 85)
mlab.show()

# histogram of neurons in the a cube
ch = ['yellow', 'c']#lepr, glp1r
fig, ax = plt.subplots()
binsj = np.int(np.amax([np.amax(density1), np.amax(scale*density2)]))+1
_,_,_ = plt.hist(density1.reshape(len(all1)*step*stepp*step), bins=binsj, color=ch[0], alpha=0.8, range=[1, binsj]); #plt.yscale('log', nonposy='clip')
_,_,_ = plt.hist(scale*density2.reshape(len(all2)*step*stepp*step), bins=binsj, color=ch[1], alpha=0.8, range=[1, binsj]); #plt.yscale('log', nonposy='clip')
plt.show()

