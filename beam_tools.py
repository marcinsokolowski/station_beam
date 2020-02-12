#!/usr/bin/python
"""Tools for transforming and plotting a primary beam generated with beam_full_EE.py.
Tim Colegate, updated from Randall Wayth's mwa_tile.py.
"""

import pdb

import numpy as np
import math
import astropy.io.fits as pyfits

import logging
logger = logging.getLogger('beam_tools')

def get_sweet_spots(sweet_spot_file='MWA_sweet_spot_gridpoints.csv'):
    """Load MWA beamformer sweet spots (including pointing direction and delays)"""
    sweet_spots={}
    temp=np.loadtxt(sweet_spot_file,delimiter=',')
    sweet_spots['astro_az']=temp[:,0]
    sweet_spots['ZA']=temp[:,1]
    sweet_spots['gridpoint']=temp[:,2]    
    sweet_spots['delays']=temp[:,3:19]        
    return sweet_spots

def gridpoint2delays(gridpoint, sweet_spot_path='MWA_sweet_spot_gridpoints.csv'):     
    """Get delays at a given gridpoint"""
    sweet_spots=get_sweet_spots(sweet_spot_path)
    gridpoint_idx=np.where(sweet_spots['gridpoint']==gridpoint)
    delays=np.vstack((sweet_spots['delays'][gridpoint_idx],
                  sweet_spots['delays'][gridpoint_idx]))
    return delays

# TEst :
# cd /home/msok/ska/aavs/aavs0.5/trunk/simulations/FEKO/beam_models/MWA_EE/MWAtools_pb/eda
# python ./beam_tools.py
#
# AZ MAP : joe out_az.txt 
# !!!!!!!!!!!!!!
# WARNING !!!! This code returns AZZA map in FEKO convention ( see out_az.txt : East=0 or 360 deg (left), North=90 deg (up), West=180 deg (right), South=270 deg (bottom) )
#     
#               North 90 deg
#  East 0 deg                   West 180 deg
#               South 270 deg
# 
# ZA MAP : joe out_za.txt
# This one looks normal with za=0 in the center and increasing symetrically as we go away from center 
def makeAZZA(npix=256,projection='SIN',azim_from_north=False,return_all=False):
    """
    Make azimuth and zenith angle arrays for a square image of side npix
    Projection is SIN or ZEA, all-sky
    Returns (az,za). Angles are in radian. ZA values are nan beyond the horizon
    """
    # build az and za arrays
    # use linspace to ensure we go to horizon on all sides
    z = np.linspace(-npix/2.0,npix/2.0,num=npix) # was dtype=np.float32  
    x = np.empty((npix,npix),dtype=np.float32)
    y = np.empty((npix,npix),dtype=np.float32)
    for i in range(npix):
        y[i,0:] = z
        x[0:,i] = z
    d = np.sqrt(x*x + y*y)/(npix/2)
    # only select pixels above horizon
    t = (d <= 1.0)
    za = np.zeros((npix,npix),dtype=np.float32)*np.NaN
    if projection=='SIN':
        za[t] = np.arcsin(d[t])
        logger.info('Using slant orthographic projection')
    elif projection=='ZEA':
        d=d*2**0.5; #ZEA requires R to extend beyond 1.
        za[t] = 2*np.arcsin(d[t]/2.0)
        logger.info('Using zenithal equal area projection')
    else:
        e='Projection %s not found'%projection
        logger.error(e)
        raise ValueError(e)
    az = np.arctan2(y,x) # to match (x,y) map and azimuth from the North it should be np.arctan2(x,y) !!!

    if azim_from_north :
        negative = (az < 0)
        az[negative] += math.pi*2.00         
    else :   
        az = az + math.pi #0 to 2pi
        az = 2*math.pi - az #Change to clockwise from top (when origin is in top-left)

    if return_all :
       return (az,za,x,y,z,d)
    else :
       return (az,za)

# WARNING : same comments as above (to makeAZZA function) apply here to a generated azza map !!!
def makeAZZA_dOMEGA(npix=256,projection='SIN'):
    """
    Make azimuth and zenith angle arrays for a square image of side npix
    Projection is SIN or ZEA, all-sky
    Returns (az,za). Angles are in radian. ZA values are nan beyond the horizon
    """
    # build az and za arrays
    # use linspace to ensure we go to horizon on all sides
    z = np.linspace(-npix/2.0,npix/2.0,num=npix) # does not work on bighorns server : ,dtype=np.float32)   
    x = np.empty((npix,npix),dtype=np.float32)
    y = np.empty((npix,npix),dtype=np.float32)
    dOMEGA = np.empty((npix,npix),dtype=np.float32)
    
    for i in range(npix):
        y[i,0:] = z
        x[0:,i] = z
    d = np.sqrt(x*x + y*y)/(npix/2)
    # only select pixels above horizon
    t = (d <= 1.0)
    n_total = t.sum()
    dOMEGA.fill( math.pi*2.00/n_total )

    za = np.zeros((npix,npix),dtype=np.float32)*np.NaN
    if projection=='SIN':
        za[t] = np.arcsin(d[t])
        dOMEGA=np.cos(za)*math.pi*2.00/n_total
        logger.info('Using slant orthographic projection')
        
    elif projection=='ZEA': # https://casa.nrao.edu/aips2_docs/memos/107/node2.html#SECTION00022200000000000000
                            # https://asdf-standard.readthedocs.io/en/stable/schemas/stsci.edu/asdf/transform/zenithal_equal_area-1.0.0.html
        d=d*2**0.5; #ZEA requires R to extend beyond 1.
        za[t] = 2*np.arcsin(d[t]/2.0) # = 2*arcsin( d * (sqrt(2)/2) ) = 2 * arcsin( d / sqrt(2) )
        logger.info('Using zenithal equal area projection')
    else:
        e='Projection %s not found'%projection
        logger.error(e)
        raise ValueError(e)
    az = np.arctan2(y,x)
    az=az+math.pi #0 to 2pi
    az=2*math.pi-az #Change to clockwise from top (when origin is in top-left)
    
    dOMEGA_sum = dOMEGA.sum()
    print "DEBUG : dOMEGA_sum = %.8f" % (dOMEGA_sum)

    return (az,za,n_total,dOMEGA)

    
def makeUnpolInstrumentalResponse(j1,j2):
    #TODO: check this description below. I think Jones dimensions are now swapped
    """
    Form the visibility matrix in instrumental response from two Jones
    matrices assuming unpolarised sources (hence the brightness matrix is
    the identity matrix)
    Input: j1,j2: Jones matrices of dimension[za][az][2][2]
    Returns: [za][az][[xx,xy],[yx,yy]] where "X" and "Y" are defined by the receptors
    of the Dipole object used in the ApertureArray. Hence to get "XX", you want
    result[za][az][0][0] and for "YY" you want result[za][az][1][1]
    """
    result = np.empty_like(j1)

    result[0,0] = j1[0,0]*j2[0,0].conjugate() + j1[0,1]*j2[0,1].conjugate()
    result[1,1] = j1[1,0]*j2[1,0].conjugate() + j1[1,1]*j2[1,1].conjugate()
    result[0,1] = j1[0,0]*j2[1,0].conjugate() + j1[0,1]*j2[1,1].conjugate()
    result[1,0] = j1[1,0]*j2[0,0].conjugate() + j1[1,1]*j2[0,1].conjugate()
    return result

def makePolInstrumentResponse(j1,j2,b):
    """
    Form the instrument response from two Jones matrices with an
    arbitrary source brightness matrix, hence arbitrary polarisation
    Returns: (xx,yy,xy,yx) where "X" and "Y" are defined by the receptors
    of the Dipole object used in the ApertureArray
    """
    # FIXME: need to work out how to do this in vectorised way.
    pass

def plotArrayJones(j, freq, filebase, title, pix_per_deg=1, j_1D=None, gridded=False):
    """
    Utility to plot the output of tile Jones matrices
    Input:
    j_1D - 1-D cut along an azimuth angle
    """
    import matplotlib.pyplot as plt
    plt.rcParams['savefig.dpi'] = 300

    for i in [0,1]:
        for ii in [0,1]:
            if j_1D is not None: #show cut
                plt.subplot(121)                
                plt.plot(np.arange(len(j_1D[i,ii]))*1.0/pix_per_deg, np.abs(j_1D[i,ii]))
                plt.title('1-D cut')
                plt.xlabel('ZA (degs)')
                plt.ylabel('magnitude')
                plt.subplot(122)

            if gridded:
                plt.imshow(np.abs(j[i,ii]), interpolation='none', extent=[0,90,360, 0])
                plt.xticks(np.arange(0, 91, 30))
                plt.yticks(np.arange(360, -1, -30))
            else:                
                plt.imshow(np.abs(j[i,ii]), interpolation='none')
            plt.suptitle('MWA %s MHz J%s%s voltage mag, %s'%(freq/1.e6,i,ii,title))
            plt.colorbar(label='magnitude')
        #    plt.gca().invert_yaxis()            
            plt.savefig('MWA_J%s%s_voltage_mag_%sMHz_%s.png'%(i,ii,freq/1.e6,filebase))
            plt.clf()

            if j_1D is not None: #show cut
                plt.subplot(121)                
                plt.plot(np.arange(len(j_1D[i,ii]))*1.0/pix_per_deg, np.angle(j_1D[i,ii])*180/math.pi)
                plt.title('1-D cut')
                plt.xlabel('ZA (deg)')
                plt.ylabel('phase (deg)')
                plt.subplot(122)

            if gridded:
                plt.imshow(np.angle(j[i,ii])*180/math.pi, interpolation='none', extent=[0,90,360, 0])
                plt.xticks(np.arange(0, 91, 30))
                plt.yticks(np.arange(360, -1, -30))            
            else:
                plt.imshow(np.angle(j[i,ii])*180/math.pi, interpolation='none')
            plt.suptitle('MWA %s MHz J%s%s voltage phase, %s'%(freq/1e6,i,ii,title))
            plt.colorbar(label='phase (deg)')
        #    plt.gca().invert_yaxis()
            plt.savefig('MWA_J%s%s_voltage_phase_%sMHz_%s.png'%(i,ii,freq/1.e6,filebase))
            plt.clf()


def exportArrayJones(j,freq,filebase):
    """
    Utility to export the output of tile Jones matrices to a .mat file
    """
    import scipy.io as io
    filename='MWA_voltage_'+str(freq/1e6)+'MHz'+filebase+'.mat'
    mydata={}

    for i in [0,1]:
        for ii in [0,1]:
            mydata['J%s%s'%(i,ii)]=j[i,ii]
    io.savemat(filename, mydata)
    
def plotVisResponse(j,freq,filebase,title,pix_per_deg,gridded=False):
    #TODO: check this description below. I think Jones dimensions are now swapped
    """
    Utility to plot the visibilty XX,YY,XY and YX response of the array for
    an unpolarised 1Jy source
    Input: j a visibility matrix (complex) of dimensions [za][az][2][2]
    """
    import matplotlib.pyplot as plt
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['axes.titlesize'] = 'medium'

    vis = makeUnpolInstrumentalResponse(j,j)
    vis=np.abs(vis)
    data={'XX':vis[0,0], 'YY':vis[1,1], 
          'XY':vis[0,1], 'YX':vis[1,0]}
    for key,val in data.iteritems():
        my_max=np.max(val,axis=1) #Max za for each az
        max_idx=np.argmax(my_max) #Find az index
        if gridded: #show cut
            plt.subplot(121)
            plt.plot(np.arange(len(j[0,0,0,:]))*1.0/pix_per_deg, val[max_idx,:])
            plt.title('Cut at Az=%.1f'%((max_idx+1)*1/pix_per_deg))
            plt.xlabel('ZA (degs)')
            plt.ylabel('power')
            plt.subplot(122)
            plt.imshow(val, interpolation='none', extent=[0,90,360, 0])
            plt.xticks(np.arange(0, 91, 30))
            plt.yticks(np.arange(360, -1, -30))
        else:
            plt.imshow(val, interpolation='none')
        plt.colorbar(label='power')
#        plt.gca().invert_yaxis()
        
        plt.suptitle('MWA '+str(freq/1.e6)+'MHz '+key+' mag '+title)    
        plt.savefig('MWA_'+key+'_mag_'+str(freq/1.e6)+'MHz_'+filebase+'.png')
        plt.clf()
        
# if floats needed :
# awk '{delays=$12;gsub(",",".,",delays);gsub("]",".]",delays);print $12" -> "delays;}' mwa_sweet_spots.txt
# awk '{delays=$12;gsub(",",".,",delays);gsub("]",".]",delays);$12=delays;print $0;}' mwa_sweet_spots.txt > mwa_sweet_spots.delays_float
#  name  , number , azimuth  , elevation ,        za        ,                     delays                     
all_grid_points = {
0  : [ 0 , 0 , 90 , 0 , [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]   ],
1  : [ 1 , 0 , 83.1912 , 6.80880000000001 , [3,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0]   ],
2  : [ 2 , 90 , 83.1912 , 6.80880000000001 , [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]   ],
3  : [ 3 , 180 , 83.1912 , 6.80880000000001 , [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]   ],
4  : [ 4 , 270 , 83.1912 , 6.80880000000001 , [3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0]   ],
5  : [ 5 , 45 , 80.348 , 9.652 , [3,4,5,6,2,3,4,5,1,2,3,4,0,1,2,3]   ],
6  : [ 6 , 135 , 80.348 , 9.652 , [0,1,2,3,1,2,3,4,2,3,4,5,3,4,5,6]   ],
7  : [ 7 , 225 , 80.348 , 9.652 , [3,2,1,0,4,3,2,1,5,4,3,2,6,5,4,3]   ],
8  : [ 8 , 315 , 80.348 , 9.652 , [6,5,4,3,5,4,3,2,4,3,2,1,3,2,1,0]   ],
9  : [ 9 , 0 , 76.2838 , 13.7162 , [6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0]   ],
10  : [ 10 , 90 , 76.2838 , 13.7162 , [0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]   ],
11  : [ 11 , 180 , 76.2838 , 13.7162 , [0,0,0,0,2,2,2,2,4,4,4,4,6,6,6,6]   ],
12  : [ 12 , 270 , 76.2838 , 13.7162 , [6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0]   ],
13  : [ 13 , 26.5651 , 74.6271 , 15.3729 , [6,7,8,9,4,5,6,7,2,3,4,5,0,1,2,3]   ],
14  : [ 14 , 63.4349 , 74.6271 , 15.3729 , [3,5,7,9,2,4,6,8,1,3,5,7,0,2,4,6]   ],
15  : [ 15 , 116.5651 , 74.6271 , 15.3729 , [0,2,4,6,1,3,5,7,2,4,6,8,3,5,7,9]   ],
16  : [ 16 , 153.4349 , 74.6271 , 15.3729 , [0,1,2,3,2,3,4,5,4,5,6,7,6,7,8,9]   ],
17  : [ 17 , 206.5651 , 74.6271 , 15.3729 , [3,2,1,0,5,4,3,2,7,6,5,4,9,8,7,6]   ],
18  : [ 18 , 243.4349 , 74.6271 , 15.3729 , [6,4,2,0,7,5,3,1,8,6,4,2,9,7,5,3]   ],
19  : [ 19 , 296.5651 , 74.6271 , 15.3729 , [9,7,5,3,8,6,4,2,7,5,3,1,6,4,2,0]   ],
20  : [ 20 , 333.4349 , 74.6271 , 15.3729 , [9,8,7,6,7,6,5,4,5,4,3,2,3,2,1,0]   ],
21  : [ 21 , 45 , 70.4075 , 19.5925 , [6,8,10,12,4,6,8,10,2,4,6,8,0,2,4,6]   ],
22  : [ 22 , 135 , 70.4075 , 19.5925 , [0,2,4,6,2,4,6,8,4,6,8,10,6,8,10,12]   ],
23  : [ 23 , 225 , 70.4075 , 19.5925 , [6,4,2,0,8,6,4,2,10,8,6,4,12,10,8,6]   ],
24  : [ 24 , 315 , 70.4075 , 19.5925 , [12,10,8,6,10,8,6,4,8,6,4,2,6,4,2,0]   ],
25  : [ 25 , 0 , 69.1655 , 20.8345 , [9,9,9,9,6,6,6,6,3,3,3,3,0,0,0,0]   ],
26  : [ 26 , 90 , 69.1655 , 20.8345 , [0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9]   ],
27  : [ 27 , 180 , 69.1655 , 20.8345 , [0,0,0,0,3,3,3,3,6,6,6,6,9,9,9,9]   ],
28  : [ 28 , 270 , 69.1655 , 20.8345 , [9,6,3,0,9,6,3,0,9,6,3,0,9,6,3,0]   ],
29  : [ 29 , 18.4349 , 67.9813 , 22.0187 , [9,10,11,12,6,7,8,9,3,4,5,6,0,1,2,3]   ],
30  : [ 30 , 71.5651 , 67.9813 , 22.0187 , [3,6,9,12,2,5,8,11,1,4,7,10,0,3,6,9]   ],
31  : [ 31 , 108.4349 , 67.9813 , 22.0187 , [0,3,6,9,1,4,7,10,2,5,8,11,3,6,9,12]   ],
32  : [ 32 , 161.5651 , 67.9813 , 22.0187 , [0,1,2,3,3,4,5,6,6,7,8,9,9,10,11,12]   ],
33  : [ 33 , 198.4349 , 67.9813 , 22.0187 , [3,2,1,0,6,5,4,3,9,8,7,6,12,11,10,9]   ],
34  : [ 34 , 251.5651 , 67.9813 , 22.0187 , [9,6,3,0,10,7,4,1,11,8,5,2,12,9,6,3]   ],
35  : [ 35 , 288.4349 , 67.9813 , 22.0187 , [12,9,6,3,11,8,5,2,10,7,4,1,9,6,3,0]   ],
36  : [ 36 , 341.5651 , 67.9813 , 22.0187 , [12,11,10,9,9,8,7,6,6,5,4,3,3,2,1,0]   ],
37  : [ 37 , 33.6901 , 64.6934 , 25.3066 , [9,11,13,15,6,8,10,12,3,5,7,9,0,2,4,6]   ],
38  : [ 38 , 56.3099 , 64.6934 , 25.3066 , [6,9,12,15,4,7,10,13,2,5,8,11,0,3,6,9]   ],
39  : [ 39 , 123.6901 , 64.6934 , 25.3066 , [0,3,6,9,2,5,8,11,4,7,10,13,6,9,12,15]   ],
40  : [ 40 , 146.3099 , 64.6934 , 25.3066 , [0,2,4,6,3,5,7,9,6,8,10,12,9,11,13,15]   ],
41  : [ 41 , 213.6901 , 64.6934 , 25.3066 , [6,4,2,0,9,7,5,3,12,10,8,6,15,13,11,9]   ],
42  : [ 42 , 236.3099 , 64.6934 , 25.3066 , [9,6,3,0,11,8,5,2,13,10,7,4,15,12,9,6]   ],
43  : [ 43 , 303.6901 , 64.6934 , 25.3066 , [15,12,9,6,13,10,7,4,11,8,5,2,9,6,3,0]   ],
44  : [ 44 , 326.3099 , 64.6934 , 25.3066 , [15,13,11,9,12,10,8,6,9,7,5,3,6,4,2,0]   ],
45  : [ 45 , 0 , 61.691 , 28.309 , [12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0]   ],
46  : [ 46 , 90 , 61.691 , 28.309 , [0,4,8,12,0,4,8,12,0,4,8,12,0,4,8,12]   ],
47  : [ 47 , 180 , 61.691 , 28.309 , [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12]   ],
48  : [ 48 , 270 , 61.691 , 28.309 , [12,8,4,0,12,8,4,0,12,8,4,0,12,8,4,0]   ],
49  : [ 49 , 14.0362 , 60.7369 , 29.2631 , [12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3]   ],
50  : [ 50 , 75.9638 , 60.7369 , 29.2631 , [3,7,11,15,2,6,10,14,1,5,9,13,0,4,8,12]   ],
51  : [ 51 , 104.0362 , 60.7369 , 29.2631 , [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]   ],
52  : [ 52 , 165.9638 , 60.7369 , 29.2631 , [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]   ],
53  : [ 53 , 194.0362 , 60.7369 , 29.2631 , [3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12]   ],
54  : [ 54 , 255.9638 , 60.7369 , 29.2631 , [12,8,4,0,13,9,5,1,14,10,6,2,15,11,7,3]   ],
55  : [ 55 , 284.0362 , 60.7369 , 29.2631 , [15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0]   ],
56  : [ 56 , 345.9638 , 60.7369 , 29.2631 , [15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]   ],
57  : [ 57 , 45 , 59.8018 , 30.1982 , [9,12,15,18,6,9,12,15,3,6,9,12,0,3,6,9]   ],
58  : [ 58 , 135 , 59.8018 , 30.1982 , [0,3,6,9,3,6,9,12,6,9,12,15,9,12,15,18]   ],
59  : [ 59 , 225 , 59.8018 , 30.1982 , [9,6,3,0,12,9,6,3,15,12,9,6,18,15,12,9]   ],
60  : [ 60 , 315 , 59.8018 , 30.1982 , [18,15,12,9,15,12,9,6,12,9,6,3,9,6,3,0]   ],
61  : [ 61 , 26.5651 , 57.981 , 32.019 , [12,14,16,18,8,10,12,14,4,6,8,10,0,2,4,6]   ],
62  : [ 62 , 63.4349 , 57.981 , 32.019 , [6,10,14,18,4,8,12,16,2,6,10,14,0,4,8,12]   ],
63  : [ 63 , 116.5651 , 57.981 , 32.019 , [0,4,8,12,2,6,10,14,4,8,12,16,6,10,14,18]   ],
64  : [ 64 , 153.4349 , 57.981 , 32.019 , [0,2,4,6,4,6,8,10,8,10,12,14,12,14,16,18]   ],
65  : [ 65 , 206.5651 , 57.981 , 32.019 , [6,4,2,0,10,8,6,4,14,12,10,8,18,16,14,12]   ],
66  : [ 66 , 243.4349 , 57.981 , 32.019 , [12,8,4,0,14,10,6,2,16,12,8,4,18,14,10,6]   ],
67  : [ 67 , 296.5651 , 57.981 , 32.019 , [18,14,10,6,16,12,8,4,14,10,6,2,12,8,4,0]   ],
68  : [ 68 , 333.4349 , 57.981 , 32.019 , [18,16,14,12,14,12,10,8,10,8,6,4,6,4,2,0]   ],
69  : [ 69 , 0 , 53.6453 , 36.3547 , [15,15,15,15,10,10,10,10,5,5,5,5,0,0,0,0]   ],
70  : [ 70 , 36.8699 , 53.6453 , 36.3547 , [12,15,18,21,8,11,14,17,4,7,10,13,0,3,6,9]   ],
71  : [ 71 , 53.1301 , 53.6453 , 36.3547 , [9,13,17,21,6,10,14,18,3,7,11,15,0,4,8,12]   ],
72  : [ 72 , 90 , 53.6453 , 36.3547 , [0,5,10,15,0,5,10,15,0,5,10,15,0,5,10,15]   ],
73  : [ 73 , 126.8699 , 53.6453 , 36.3547 , [0,4,8,12,3,7,11,15,6,10,14,18,9,13,17,21]   ],
74  : [ 74 , 143.1301 , 53.6453 , 36.3547 , [0,3,6,9,4,7,10,13,8,11,14,17,12,15,18,21]   ],
75  : [ 75 , 180 , 53.6453 , 36.3547 , [0,0,0,0,5,5,5,5,10,10,10,10,15,15,15,15]   ],
76  : [ 76 , 216.8699 , 53.6453 , 36.3547 , [9,6,3,0,13,10,7,4,17,14,11,8,21,18,15,12]   ],
77  : [ 77 , 233.1301 , 53.6453 , 36.3547 , [12,8,4,0,15,11,7,3,18,14,10,6,21,17,13,9]   ],
78  : [ 78 , 270 , 53.6453 , 36.3547 , [15,10,5,0,15,10,5,0,15,10,5,0,15,10,5,0]   ],
79  : [ 79 , 306.8699 , 53.6453 , 36.3547 , [21,17,13,9,18,14,10,6,15,11,7,3,12,8,4,0]   ],
80  : [ 80 , 323.1301 , 53.6453 , 36.3547 , [21,18,15,12,17,14,11,8,13,10,7,4,9,6,3,0]   ],
81  : [ 81 , 11.3099 , 52.8056 , 37.1944 , [15,16,17,18,10,11,12,13,5,6,7,8,0,1,2,3]   ],
82  : [ 82 , 78.6901 , 52.8056 , 37.1944 , [3,8,13,18,2,7,12,17,1,6,11,16,0,5,10,15]   ],
83  : [ 83 , 101.3099 , 52.8056 , 37.1944 , [0,5,10,15,1,6,11,16,2,7,12,17,3,8,13,18]   ],
84  : [ 84 , 168.6901 , 52.8056 , 37.1944 , [0,1,2,3,5,6,7,8,10,11,12,13,15,16,17,18]   ],
85  : [ 85 , 191.3099 , 52.8056 , 37.1944 , [3,2,1,0,8,7,6,5,13,12,11,10,18,17,16,15]   ],
86  : [ 86 , 258.6901 , 52.8056 , 37.1944 , [15,10,5,0,16,11,6,1,17,12,7,2,18,13,8,3]   ],
87  : [ 87 , 281.3099 , 52.8056 , 37.1944 , [18,13,8,3,17,12,7,2,16,11,6,1,15,10,5,0]   ],
88  : [ 88 , 348.6901 , 52.8056 , 37.1944 , [18,17,16,15,13,12,11,10,8,7,6,5,3,2,1,0]   ],
89  : [ 89 , 21.8014 , 50.3239 , 39.6761 , [15,17,19,21,10,12,14,16,5,7,9,11,0,2,4,6]   ],
90  : [ 90 , 68.1986 , 50.3239 , 39.6761 , [6,11,16,21,4,9,14,19,2,7,12,17,0,5,10,15]   ],
91  : [ 91 , 111.8014 , 50.3239 , 39.6761 , [0,5,10,15,2,7,12,17,4,9,14,19,6,11,16,21]   ],
92  : [ 92 , 158.1986 , 50.3239 , 39.6761 , [0,2,4,6,5,7,9,11,10,12,14,16,15,17,19,21]   ],
93  : [ 93 , 201.8014 , 50.3239 , 39.6761 , [6,4,2,0,11,9,7,5,16,14,12,10,21,19,17,15]   ],
94  : [ 94 , 248.1986 , 50.3239 , 39.6761 , [15,10,5,0,17,12,7,2,19,14,9,4,21,16,11,6]   ],
95  : [ 95 , 291.8014 , 50.3239 , 39.6761 , [21,16,11,6,19,14,9,4,17,12,7,2,15,10,5,0]   ],
96  : [ 96 , 338.1986 , 50.3239 , 39.6761 , [21,19,17,15,16,14,12,10,11,9,7,5,6,4,2,0]   ],
97  : [ 97 , 45 , 47.8822 , 42.1178 , [12,16,20,24,8,12,16,20,4,8,12,16,0,4,8,12]   ],
98  : [ 98 , 135 , 47.8822 , 42.1178 , [0,4,8,12,4,8,12,16,8,12,16,20,12,16,20,24]   ],
99  : [ 99 , 225 , 47.8822 , 42.1178 , [12,8,4,0,16,12,8,4,20,16,12,8,24,20,16,12]   ],
100  : [ 100 , 315 , 47.8822 , 42.1178 , [24,20,16,12,20,16,12,8,16,12,8,4,12,8,4,0]   ],
101  : [ 101 , 30.9638 , 46.2671 , 43.7329 , [15,18,21,24,10,13,16,19,5,8,11,14,0,3,6,9]   ],
102  : [ 102 , 59.0362 , 46.2671 , 43.7329 , [9,14,19,24,6,11,16,21,3,8,13,18,0,5,10,15]   ],
103  : [ 103 , 120.9638 , 46.2671 , 43.7329 , [0,5,10,15,3,8,13,18,6,11,16,21,9,14,19,24]   ],
104  : [ 104 , 149.0362 , 46.2671 , 43.7329 , [0,3,6,9,5,8,11,14,10,13,16,19,15,18,21,24]   ],
105  : [ 105 , 210.9638 , 46.2671 , 43.7329 , [9,6,3,0,14,11,8,5,19,16,13,10,24,21,18,15]   ],
106  : [ 106 , 239.0362 , 46.2671 , 43.7329 , [15,10,5,0,18,13,8,3,21,16,11,6,24,19,14,9]   ],
107  : [ 107 , 300.9638 , 46.2671 , 43.7329 , [24,19,14,9,21,16,11,6,18,13,8,3,15,10,5,0]   ],
108  : [ 108 , 329.0362 , 46.2671 , 43.7329 , [24,21,18,15,19,16,13,10,14,11,8,5,9,6,3,0]   ],
109  : [ 109 , 0 , 44.656 , 45.344 , [18,18,18,18,12,12,12,12,6,6,6,6,0,0,0,0]   ],
110  : [ 110 , 90 , 44.656 , 45.344 , [0,6,12,18,0,6,12,18,0,6,12,18,0,6,12,18]   ],
111  : [ 111 , 180 , 44.656 , 45.344 , [0,0,0,0,6,6,6,6,12,12,12,12,18,18,18,18]   ],
112  : [ 112 , 270 , 44.656 , 45.344 , [18,12,6,0,18,12,6,0,18,12,6,0,18,12,6,0]   ],
113  : [ 113 , 9.4623 , 43.8504 , 46.1496 , [18,19,20,21,12,13,14,15,6,7,8,9,0,1,2,3]   ],
114  : [ 114 , 80.5377 , 43.8504 , 46.1496 , [3,9,15,21,2,8,14,20,1,7,13,19,0,6,12,18]   ],
115  : [ 115 , 99.4623 , 43.8504 , 46.1496 , [0,6,12,18,1,7,13,19,2,8,14,20,3,9,15,21]   ],
116  : [ 116 , 170.5377 , 43.8504 , 46.1496 , [0,1,2,3,6,7,8,9,12,13,14,15,18,19,20,21]   ],
117  : [ 117 , 189.4623 , 43.8504 , 46.1496 , [3,2,1,0,9,8,7,6,15,14,13,12,21,20,19,18]   ],
118  : [ 118 , 260.5377 , 43.8504 , 46.1496 , [18,12,6,0,19,13,7,1,20,14,8,2,21,15,9,3]   ],
119  : [ 119 , 279.4623 , 43.8504 , 46.1496 , [21,15,9,3,20,14,8,2,19,13,7,1,18,12,6,0]   ],
120  : [ 120 , 350.5377 , 43.8504 , 46.1496 , [21,20,19,18,15,14,13,12,9,8,7,6,3,2,1,0]   ],
121  : [ 121 , 18.4349 , 41.4255 , 48.5745 , [18,20,22,24,12,14,16,18,6,8,10,12,0,2,4,6]   ],
122  : [ 122 , 71.5651 , 41.4255 , 48.5745 , [6,12,18,24,4,10,16,22,2,8,14,20,0,6,12,18]   ],
123  : [ 123 , 108.4349 , 41.4255 , 48.5745 , [0,6,12,18,2,8,14,20,4,10,16,22,6,12,18,24]   ],
124  : [ 124 , 161.5651 , 41.4255 , 48.5745 , [0,2,4,6,6,8,10,12,12,14,16,18,18,20,22,24]   ],
125  : [ 125 , 198.4349 , 41.4255 , 48.5745 , [6,4,2,0,12,10,8,6,18,16,14,12,24,22,20,18]   ],
126  : [ 126 , 251.5651 , 41.4255 , 48.5745 , [18,12,6,0,20,14,8,2,22,16,10,4,24,18,12,6]   ],
127  : [ 127 , 288.4349 , 41.4255 , 48.5745 , [24,18,12,6,22,16,10,4,20,14,8,2,18,12,6,0]   ],
128  : [ 128 , 341.5651 , 41.4255 , 48.5745 , [24,22,20,18,18,16,14,12,12,10,8,6,6,4,2,0]   ],
129  : [ 129 , 38.6598 , 40.6123 , 49.3877 , [15,19,23,27,10,14,18,22,5,9,13,17,0,4,8,12]   ],
130  : [ 130 , 51.3402 , 40.6123 , 49.3877 , [12,17,22,27,8,13,18,23,4,9,14,19,0,5,10,15]   ],
131  : [ 131 , 128.6598 , 40.6123 , 49.3877 , [0,5,10,15,4,9,14,19,8,13,18,23,12,17,22,27]   ],
132  : [ 132 , 141.3402 , 40.6123 , 49.3877 , [0,4,8,12,5,9,13,17,10,14,18,22,15,19,23,27]   ],
133  : [ 133 , 218.6598 , 40.6123 , 49.3877 , [12,8,4,0,17,13,9,5,22,18,14,10,27,23,19,15]   ],
134  : [ 134 , 231.3402 , 40.6123 , 49.3877 , [15,10,5,0,19,14,9,4,23,18,13,8,27,22,17,12]   ],
135  : [ 135 , 308.6598 , 40.6123 , 49.3877 , [27,22,17,12,23,18,13,8,19,14,9,4,15,10,5,0]   ],
136  : [ 136 , 321.3402 , 40.6123 , 49.3877 , [27,23,19,15,22,18,14,10,17,13,9,5,12,8,4,0]   ],
137  : [ 137 , 26.5651 , 37.3163 , 52.6837 , [18,21,24,27,12,15,18,21,6,9,12,15,0,3,6,9]   ],
138  : [ 138 , 63.4349 , 37.3163 , 52.6837 , [9,15,21,27,6,12,18,24,3,9,15,21,0,6,12,18]   ],
139  : [ 139 , 116.5651 , 37.3163 , 52.6837 , [0,6,12,18,3,9,15,21,6,12,18,24,9,15,21,27]   ],
140  : [ 140 , 153.4349 , 37.3163 , 52.6837 , [0,3,6,9,6,9,12,15,12,15,18,21,18,21,24,27]   ],
141  : [ 141 , 206.5651 , 37.3163 , 52.6837 , [9,6,3,0,15,12,9,6,21,18,15,12,27,24,21,18]   ],
142  : [ 142 , 243.4349 , 37.3163 , 52.6837 , [18,12,6,0,21,15,9,3,24,18,12,6,27,21,15,9]   ],
143  : [ 143 , 296.5651 , 37.3163 , 52.6837 , [27,21,15,9,24,18,12,6,21,15,9,3,18,12,6,0]   ],
144  : [ 144 , 333.4349 , 37.3163 , 52.6837 , [27,24,21,18,21,18,15,12,15,12,9,6,9,6,3,0]   ],
145  : [ 145 , 0 , 33.912 , 56.088 , [21,21,21,21,14,14,14,14,7,7,7,7,0,0,0,0]   ],
146  : [ 146 , 90 , 33.912 , 56.088 , [0,7,14,21,0,7,14,21,0,7,14,21,0,7,14,21]   ],
147  : [ 147 , 180 , 33.912 , 56.088 , [0,0,0,0,7,7,7,7,14,14,14,14,21,21,21,21]   ],
148  : [ 148 , 270 , 33.912 , 56.088 , [21,14,7,0,21,14,7,0,21,14,7,0,21,14,7,0]   ],
149  : [ 149 , 8.1301 , 33.0368 , 56.9632 , [21,22,23,24,14,15,16,17,7,8,9,10,0,1,2,3]   ],
150  : [ 150 , 45 , 33.0368 , 56.9632 , [15,20,25,30,10,15,20,25,5,10,15,20,0,5,10,15]   ],
151  : [ 151 , 81.8699 , 33.0368 , 56.9632 , [3,10,17,24,2,9,16,23,1,8,15,22,0,7,14,21]   ],
152  : [ 152 , 98.1301 , 33.0368 , 56.9632 , [0,7,14,21,1,8,15,22,2,9,16,23,3,10,17,24]   ],
153  : [ 153 , 135 , 33.0368 , 56.9632 , [0,5,10,15,5,10,15,20,10,15,20,25,15,20,25,30]   ],
154  : [ 154 , 171.8699 , 33.0368 , 56.9632 , [0,1,2,3,7,8,9,10,14,15,16,17,21,22,23,24]   ],
155  : [ 155 , 188.1301 , 33.0368 , 56.9632 , [3,2,1,0,10,9,8,7,17,16,15,14,24,23,22,21]   ],
156  : [ 156 , 225 , 33.0368 , 56.9632 , [15,10,5,0,20,15,10,5,25,20,15,10,30,25,20,15]   ],
157  : [ 157 , 261.8699 , 33.0368 , 56.9632 , [21,14,7,0,22,15,8,1,23,16,9,2,24,17,10,3]   ],
158  : [ 158 , 278.1301 , 33.0368 , 56.9632 , [24,17,10,3,23,16,9,2,22,15,8,1,21,14,7,0]   ],
159  : [ 159 , 315 , 33.0368 , 56.9632 , [30,25,20,15,25,20,15,10,20,15,10,5,15,10,5,0]   ],
160  : [ 160 , 351.8699 , 33.0368 , 56.9632 , [24,23,22,21,17,16,15,14,10,9,8,7,3,2,1,0]   ],
161  : [ 161 , 33.6901 , 31.2488 , 58.7512 , [18,22,26,30,12,16,20,24,6,10,14,18,0,4,8,12]   ],
162  : [ 162 , 56.3099 , 31.2488 , 58.7512 , [12,18,24,30,8,14,20,26,4,10,16,22,0,6,12,18]   ],
163  : [ 163 , 123.6901 , 31.2488 , 58.7512 , [0,6,12,18,4,10,16,22,8,14,20,26,12,18,24,30]   ],
164  : [ 164 , 146.3099 , 31.2488 , 58.7512 , [0,4,8,12,6,10,14,18,12,16,20,24,18,22,26,30]   ],
165  : [ 165 , 213.6901 , 31.2488 , 58.7512 , [12,8,4,0,18,14,10,6,24,20,16,12,30,26,22,18]   ],
166  : [ 166 , 236.3099 , 31.2488 , 58.7512 , [18,12,6,0,22,16,10,4,26,20,14,8,30,24,18,12]   ],
167  : [ 167 , 303.6901 , 31.2488 , 58.7512 , [30,24,18,12,26,20,14,8,22,16,10,4,18,12,6,0]   ],
168  : [ 168 , 326.3099 , 31.2488 , 58.7512 , [30,26,22,18,24,20,16,12,18,14,10,6,12,8,4,0]   ],
169  : [ 169 , 15.9454 , 30.3331 , 59.6669 , [21,23,25,27,14,16,18,20,7,9,11,13,0,2,4,6]   ],
170  : [ 170 , 74.0546 , 30.3331 , 59.6669 , [6,13,20,27,4,11,18,25,2,9,16,23,0,7,14,21]   ],
171  : [ 171 , 105.9454 , 30.3331 , 59.6669 , [0,7,14,21,2,9,16,23,4,11,18,25,6,13,20,27]   ],
172  : [ 172 , 164.0546 , 30.3331 , 59.6669 , [0,2,4,6,7,9,11,13,14,16,18,20,21,23,25,27]   ],
173  : [ 173 , 195.9454 , 30.3331 , 59.6669 , [6,4,2,0,13,11,9,7,20,18,16,14,27,25,23,21]   ],
174  : [ 174 , 254.0546 , 30.3331 , 59.6669 , [21,14,7,0,23,16,9,2,25,18,11,4,27,20,13,6]   ],
175  : [ 175 , 285.9454 , 30.3331 , 59.6669 , [27,20,13,6,25,18,11,4,23,16,9,2,21,14,7,0]   ],
176  : [ 176 , 344.0546 , 30.3331 , 59.6669 , [27,25,23,21,20,18,16,14,13,11,9,7,6,4,2,0]   ],
177  : [ 177 , 23.1986 , 25.4582 , 64.5418 , [21,24,27,30,14,17,20,23,7,10,13,16,0,3,6,9]   ],
178  : [ 178 , 66.8014 , 25.4582 , 64.5418 , [9,16,23,30,6,13,20,27,3,10,17,24,0,7,14,21]   ],
179  : [ 179 , 113.1986 , 25.4582 , 64.5418 , [0,7,14,21,3,10,17,24,6,13,20,27,9,16,23,30]   ],
180  : [ 180 , 156.8014 , 25.4582 , 64.5418 , [0,3,6,9,7,10,13,16,14,17,20,23,21,24,27,30]   ],
181  : [ 181 , 203.1986 , 25.4582 , 64.5418 , [9,6,3,0,16,13,10,7,23,20,17,14,30,27,24,21]   ],
182  : [ 182 , 246.8014 , 25.4582 , 64.5418 , [21,14,7,0,24,17,10,3,27,20,13,6,30,23,16,9]   ],
183  : [ 183 , 293.1986 , 25.4582 , 64.5418 , [30,23,16,9,27,20,13,6,24,17,10,3,21,14,7,0]   ],
184  : [ 184 , 336.8014 , 25.4582 , 64.5418 , [30,27,24,21,23,20,17,14,16,13,10,7,9,6,3,0]   ],
185  : [ 185 , 0 , 18.4768 , 71.5232 , [24,24,24,24,16,16,16,16,8,8,8,8,0,0,0,0]   ],
186  : [ 186 , 90 , 18.4768 , 71.5232 , [0,8,16,24,0,8,16,24,0,8,16,24,0,8,16,24]   ],
187  : [ 187 , 180 , 18.4768 , 71.5232 , [0,0,0,0,8,8,8,8,16,16,16,16,24,24,24,24]   ],
188  : [ 188 , 270 , 18.4768 , 71.5232 , [24,16,8,0,24,16,8,0,24,16,8,0,24,16,8,0]   ],
189  : [ 189 , 7.125 , 17.0922 , 72.9078 , [24,25,26,27,16,17,18,19,8,9,10,11,0,1,2,3]   ],
190  : [ 190 , 82.875 , 17.0922 , 72.9078 , [3,11,19,27,2,10,18,26,1,9,17,25,0,8,16,24]   ],
191  : [ 191 , 97.125 , 17.0922 , 72.9078 , [0,8,16,24,1,9,17,25,2,10,18,26,3,11,19,27]   ],
192  : [ 192 , 172.875 , 17.0922 , 72.9078 , [0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27]   ],
193  : [ 193 , 187.125 , 17.0922 , 72.9078 , [3,2,1,0,11,10,9,8,19,18,17,16,27,26,25,24]   ],
194  : [ 194 , 262.875 , 17.0922 , 72.9078 , [24,16,8,0,25,17,9,1,26,18,10,2,27,19,11,3]   ],
195  : [ 195 , 277.125 , 17.0922 , 72.9078 , [27,19,11,3,26,18,10,2,25,17,9,1,24,16,8,0]   ],
196  : [ 196 , 352.875 , 17.0922 , 72.9078 , [27,26,25,24,19,18,17,16,11,10,9,8,3,2,1,0]   ]
}


def find_closest_gridpoint(az_deg,za_deg,gridpoint_list=all_grid_points) :
   best=gridpoint_list[0]
   deg2rad=math.pi/180.00;
   
   min_dist=1e6
   for idx in gridpoint_list :
      gridpoint=gridpoint_list[idx]
      az_gridpoint=gridpoint[1]
      alt_gridpoint=gridpoint[2]
      za_gridpoint=gridpoint[3]
      
      alt_deg=90.00-za_deg
      diff_az=az_gridpoint-az_deg

#      dist=math.fabs(za-za_gridpoint)
      cos_value = math.sin(alt_gridpoint*deg2rad)*math.sin(alt_deg*deg2rad) + math.cos(alt_gridpoint*deg2rad)*math.cos(alt_deg*deg2rad)*math.cos( diff_az*deg2rad );
      dist=math.acos(cos_value)
      if dist < min_dist :
         min_dist = dist
         best = gridpoint
         
   return best

def get_delays(gridpoint_idx,gridpoint_list=all_grid_points) :

   for idx in gridpoint_list:
      gridpoint=gridpoint_list[idx]
      if gridpoint[0] == gridpoint_idx :
         delays = gridpoint[4]
         return np.vstack((delays,delays))
      
   return None

if __name__ == "__main__":
   npix=64
   (az,za) = makeAZZA( npix , 'SIN' , azim_from_north=True )           
   
   out_f_za = open("out_za.txt", "w" )
   out_f_az = open("out_az.txt", "w" )
   
   for y in range(npix-1,-1,-1):
      line_az = ""
      line_za = ""
      for x in range(0,npix) :
         az_str = "%03.1f " % (az[x,y]*(180.00/math.pi))
         za_str = "%03.1f " % (za[x,y]*(180.00/math.pi))
         line_az += az_str
         line_za += za_str
         
      out_f_az.write( line_az + "\n" )
      out_f_za.write( line_za + "\n" )
      
         
           
   out_f_za.close()
   out_f_az.close()
 
        
   hdu = pyfits.PrimaryHDU()
   hdu.data = az*(180.0/math.pi)
   out_fits_name = "test_azza_AZ_SIN.fits"
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_fits_name , clobber=True )

   hdu.data = za*(180.0/math.pi)
   out_fits_name = "test_azza_ZA_SIN.fits"
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_fits_name , clobber=True )
   
   
   (az,za) = makeAZZA( npix , 'ZEA' )
   hdu = pyfits.PrimaryHDU()
   hdu.data = az*(180.0/math.pi)
   out_fits_name = "test_azza_AZ_ZEA.fits"
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_fits_name , clobber=True )

   hdu.data = za*(180.0/math.pi)
   out_fits_name = "test_azza_ZA_ZEA.fits"
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_fits_name , clobber=True )

   
      
         