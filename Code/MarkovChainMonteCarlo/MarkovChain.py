# This code generates stepping locations using a Markov Chain Monte Carlo model.
# The output of the code are histograms of slope of the stepping model under two
# stepping schemes; directed scheme where the runner aims to land on flatter patches 
# of the terrain, and the uniformly random blind scheme. These are compared to experimental
# data on stepping locations that I collected from human subject experiments.


import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
import random as rnd
import math
from scipy import stats

def noise_gen(xland,yland,Dely,Delx):
	# generating von Mises noise around landing location

	kappa = 1

	xout = (rnd.vonmisesvariate(math.pi,kappa)-math.pi)*(Delx/(2*math.pi))+xland
	yout = (rnd.vonmisesvariate(math.pi,kappa)-math.pi)*(Dely/(2*math.pi))+yland

	return [xout, yout]

def footplace(xland,yland,terrht,terrstd,xed,yed,sl,sw,xprev):
	
	#define size of search space
	delx = 2*sw
	dely = 2*sl

	#define bounds of search space centered around landing spot
	if (yland - dely > yed[0]):
		ylow = max(np.argwhere(yed<yland-dely))
	else:
		ylow = np.array([0])
	if (yland + dely < yed[-1]):
		ymax = max(np.argwhere(yed<yland+dely))
	else:
		ymax = np.array([len(yed)-1])

	# move left or right depending on previous step
	if (xland > xprev):
		#print(xprev, xed[0])
		if (xland - delx > xed[0]):
			x1 = max(np.argwhere(xed<=xland-delx))
		else:
			x1 = 0
		if(xprev > xed[0]):
			x2 = max(np.argwhere(xed<=xprev))
		else:
			x2 = 0
		xlow = np.array([max(x1,x2)])
		if (xed[-1] > xland + delx):
			xmax = min(np.argwhere(xed > xland+delx))
			xmax = min(xmax,np.array([len(xed)-1]))
		else:
			xmax = np.array([len(xed)-1])
	else:
		temp = np.argwhere(xed>=xprev)
		if (not any(temp)):
			temp = np.array([len(xed)])
		xmax = min(temp[0],np.array([len(xed)-1]))
		if (any(np.argwhere(xed<xland-delx))):
			xlow = max(np.argwhere(xed<xland-delx))
		else:
			xlow = np.array([0])

	# find location of minimum slope patch
	xlow = np.asscalar(xlow)
	xmax = np.asscalar(xmax)
	ylow = np.asscalar(ylow)
	ymax = np.asscalar(ymax)
	indy = np.argmin(np.min(terrstd[0,xlow:xmax+1,ylow:ymax+1],axis=0))
	indx = np.argmin(np.min(terrstd[0,xlow:xmax+1,ylow:ymax+1],axis=1))
	
	# set stepping location to minimum patch
	if (indx + xlow < len(xed)-1):
		xland1 = (xed[indx+xlow]+xed[indx+xlow+1])/2
	else:
		xland1 = (xed[-2]+xed[-1])/2
	if (xland < xprev and xland1 < xprev):
		xland = xland1
	elif(xland < xprev and xland1 > xprev):
		xland = max(xprev - sw,xed[0])
	elif(xland > xprev and xland1 > xprev):
		xland = xland1
	else:
		xland = min(xprev + sw,xed[-1])
	if (indy+ylow < len(yed)-1):
		yland = (yed[indy+ylow]+yed[indy+ylow+1])/2
	else:
		yland = (yed[indy+ylow]+yed[-1])/2

	# add noise to stepping location
	fnout = noise_gen(xland,yland,sl,sw)
	xland = fnout[0]
	yland = fnout[1]
	if (xland < xed[0]):
		xland = xed[0]
	elif(xland>xed[-1]):
		xland = xed[-1]

	# find height and slope of terrain at landing location
	if(yland<yed[-1] and yland>yed[0] and xland > xed[0] and xland < xed[-1]):
		xind = np.max(np.argwhere(xed<xland))
		yind = np.max(np.argwhere(yed<yland))

		xind = np.asscalar(xind)
		yind = np.asscalar(yind)
		vh = terrht[0,xind,yind]
		vstd = terrstd[0,xind,yind]
	else:
		vh = []
		vstd = []

	return [vh, vstd, xland, yland]


# Importing data from .mat file

mat_contents = spio.loadmat('grid-uneven.mat')
# Arrays
htall = mat_contents['htallm'] # 10 x 21 x 23
slopeall = mat_contents['slopeallm'] # 10 x 21 x 23
tbnds = mat_contents['tbnds'] # 10 x 2
# Scalars
xhi = mat_contents['xhi']
xlow = mat_contents['xlow']

# Defining stepping parameters

meansl = 1129
stdsl = 62
meansw = 53
swstd = 18
heelx = 30
heely = 80

# Define track and simulation parameters

trackl = 25000
trackw = xhi-xlow
maxsteps = 1000

# initialize arrays

htdist = np.array([])
stdhtdist = np.array([])
bigsl =  np.zeros(maxsteps-1)
bigsw = np.zeros(maxsteps-1)
allht = np.array([])
stdht = np.array([])
x = np.zeros(maxsteps)
y = np.zeros(maxsteps)

# Initial conditions

xmid = (xlow+xhi)/2
x0 = xmid+meansw/2 
y0 = trackl/2
k = 1
x[0] = x0
y[0] = y0

# Find ends of the track
ymin = min(tbnds[:,0])
ymax = max(tbnds[:,1])

for i in range(1,maxsteps):

	# if reach beginning or end of the track, turn around and reset x position
	if (y[i-1]+((-1)**k)*meansl < -trackl/2):
		y[i] = -trackl/2
		x[i] = xmid + ((-1)**i)*meansw/2
		k = k+1
	elif (y[i-1]+((-1)**k)*meansl > trackl/2):
		y[i] = trackl/2
		x[i] = xmid + ((-1)**i)*meansw/2
		k = k+1
	else:
	# take an average step forward
		x[i] = x[i-1]+((-1)**i)*meansw
		if(x[i]<xlow):
			x[i] = xlow
		elif(x[i]>xhi):
			x[i] = xhi
		y[i] = y[i-1]+((-1)**k)*meansl

		# find which panel you are in
		if (y[i] > ymin and y[i] < ymax):
			ypos = np.argwhere(tbnds[:,1]>y[i])
			ypos = ypos[-1]

	    	#if step within capture volume
			if (tbnds[ypos,0] < y[i]):

				# load terrain properties of panel
				terrainhts = htall[ypos,:,:]
				terrainstd = slopeall[ypos,:,:]

				#create heel size grid
				xed = np.arange(xlow,xhi,heelx)
				yed = np.arange(tbnds[ypos,0],tbnds[ypos,1],heely)

				# find optimal stepping position
				returns = footplace(x[i],y[i],terrainhts,terrainstd,xed,yed,stdsl,swstd,x[i-1])
				ht = returns[0]
				hstd = returns[1]
				x[i] = returns[2]
				y[i] = returns[3]

				# save terrain properties of stepping location
				htdist = np.append(htdist,ht)
				stdhtdist = np.append(stdhtdist,hstd)

			#if step outside capture volume
			else:
				# add noise to stepping location
				newpos = noise_gen(x[i],y[i],stdsl,swstd)
				x[i] = newpos[0]
				if (x[i] < xlow):
					x[i] = xlow
				elif(x[i]>xhi):
					x[i] = xhi
				y[i] = newpos[1]
		else:
			newpos = noise_gen(x[i],y[i],stdsl,swstd)
			x[i] = newpos[0]
			if (x[i] < xlow):
				x[i] = xlow
			elif(x[i]>xhi):
				x[i] = xhi
			y[i] = newpos[1]

	   #store step length and step width     
		bigsl[i-1] = abs(y[i]-y[i-1])
		bigsw[i-1] = abs(x[i]-x[i-1])

#initialize blind runner terrain properties
htterr = np.array([])
stdhtterr = np.array([])

# store blind runner terrain properties in list
for i in range(10):
	htterr = np.append(htterr,htall[i,:,:].flatten())
	stdhtterr = np.append(stdhtterr,slopeall[i,:,:].flatten())

#Remove nans
htdist = htdist[~np.isnan(htdist)]
stdhtdist = stdhtdist[~np.isnan(stdhtdist)]
stdhtdist[stdhtdist.astype(bool)]
htterr = htterr[~np.isnan(htterr)]
stdhtterr = stdhtterr[~np.isnan(stdhtterr)]

#load regression analysis data
regressdata = spio.loadmat('regressionstuff.mat')
problanding = regressdata['landings']
terrainrough = regressdata['terrstdland']

#cleaning data
allnans = np.isnan(terrainrough)
terrainrough = terrainrough[~allnans]
problanding = problanding[~allnans]
problanding = np.true_divide(problanding,40)

#performing linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(terrainrough,problanding)

#plotting regression
xt = np.linspace(min(terrainrough),max(terrainrough),100)
line = np.ones(len(xt))*intercept + slope*xt
plt.figure(0)
plt.scatter(terrainrough,problanding)
plt.plot(xt,line,'r',label='y={:2f}x+{:2f}'.format(slope,intercept))
plt.legend()
plt.xlabel('roughness')
plt.ylabel('probability of stepping')
plt.show(0)

# load experimental data
exp_contents = spio.loadmat('footplace-un1-all.mat')
exp_slope = exp_contents['slun1']
exp_slope = exp_slope[0,:]
exp_slope = exp_slope[~np.isnan(exp_slope)]

# Plot histograms
plt.figure(1)
n1, bins1, patches1 = plt.hist(stdhtdist,25,density = True,facecolor='b',alpha=0.4,label = 'directed scheme')
n2, bin2, patches2 = plt.hist(stdhtterr,bins=bins1,density = True,facecolor='y',alpha=0.4,label = 'blind scheme')
n3, bin3, patches3 = plt.hist(exp_slope,bins=bins1,density = True,facecolor='r',alpha=0.4,label= 'experimental data')
plt.xlabel('roughness')
plt.ylabel('pdf')
plt.legend()
plt.show(1)