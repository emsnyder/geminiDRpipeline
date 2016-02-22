import numpy as np
import pylab
import pyfits
import sys
import glob
import os
from pyraf import iraf
import glob
from astropy.io import fits as fits

print '----------------------------------------------'
print 'Welcome to the Gemini data reduction pipeline!'
print '----------------------------------------------'
print ' '
print ' '

#find where user is working

mainpath = os.getcwd()

foldername=raw_input('What is the folder you are working in?    ')

# collect all fits files in the working directory

files = glob.glob('S*.fits')

sz = np.size(files)

wave     = np.chararray(sz,itemsize=10)
obj      = np.chararray(sz,itemsize=10)
grat     = np.chararray(sz,itemsize=10)
mask     = np.chararray(sz,itemsize=10)
obstype  = np.chararray(sz,itemsize=10)
obsclass = np.chararray(sz,itemsize=10)
date     = np.chararray(sz,itemsize=10)

for i in np.arange(0,sz):
    tmp         = pyfits.open(files[i])
    hdr         = tmp[0].header
    wave[i]     = (hdr['centwave'])
    obj[i]      = (hdr['object'])
    grat[i]     = (hdr['grating'])
    mask[i]     = (hdr['maskname'])
    obstype[i]  = (hdr['obstype'])
    obsclass[i] = (hdr['obsclass'])
    date[i]     = (hdr['date'])

if mask[0] == 'IFU-R': setup = 'blue'
if mask[0] == 'IFU-2': setup = 'red'

for i in np.arange(0,sz):
    print files[i],obj[i],grat[i],mask[i],obstype[i],obsclass[i],wave[i]

#select either the red or blue setups

if setup == 'blue':
    waves = np.array(['480.0', '490.0'])
    slit = 'red'
if setup == 'red':
    waves = np.array(['635.0', '640.0'])
    slit = 'both'

#pull out acquisition images

selacq = (obsclass == 'acq')
acqimages = np.extract(selacq, files)

#pull out science images

selsci1 = (obsclass == 'science') & (wave == waves[0])
sciimage1 = np.extract(selsci1,files)
n = np.size(sciimage1)
if n <= 1:
    sci11 = sciimage1[0]
else:
    sci11 = sciimage1[0]
    sci12 = sciimage1[1]
selsci2 = (obsclass == 'science') & (wave == waves[1])
sciimage2 = np.extract(selsci2,files)
if n <= 1:
    sci21 = sciimage2[0]
else:
    sci21 = sciimage2[0]
    sci22 = sciimage2[1]

#pull out arcs

selarc1 = (wave == waves[0]) & (obstype == 'ARC')
arc1 = np.extract(selarc1, files)
selarc2 = (wave == waves[1]) & (obstype == 'ARC')
arc2 = np.extract(selarc2, files)

#pull out flats

selflat1 = (wave == waves[0]) & (obstype == 'FLAT')
flat1 = np.extract(selflat1, files)[0]
flat1base = flat1.split('.')[0]
selflat2 = (wave == waves[1]) & (obstype == 'FLAT')
flat2 = np.extract(selflat2, files)[0]
flat2base = flat2.split('.')[0]

#pull out date

obsdate = date[0] #date should be same for all files

if obsdate == '2015-11-08': biasim='20151108bias.fits'

#pull out bias

#pull out standard star info

#if date = 2015 and setup = blue then standard star folder = ??
 

#pull out galaxy name

selname = (obsclass == 'science') & (wave == waves[0])
galname = np.extract(selname, obj)

#f = open('reduction.cl', 'w')

#defining all reduction tasks

def run_gfreduce(inimage,slit,cr,wt,ss,inter,os,tr,bi,):
    iraf.gemini.gmos.gfreduce(inimage, fl_gscrrej=cr, fl_wavtran=wt, fl_skysub=ss, fl_inter=inter, slits=slit, fl_over=yes, rawpath=foldername, logfile=galname[0]+'.log', bias=biasim)

def run_arcs(inimage,flatname):
    iraf.gemini.gmos.gfreduce(inimage, wavtran=no, fl_inter=no, ref='erg'+flatname, recenter=no, trace=no, fl_skysub=no, fl_gscrrej=no, fl_bias=no, fl_over=no, order=1, weights=none, slits=slit)

def plot_spaxels(inimage,verinfo):
    iraf.gemini.gmos.gfdisplay(inimage, ver=verinfo)

def wave_cal(inimage):
    iraf.gemini.gmos.gswavelength(inimage, fl_inter=yes, nlost=10)
#add coordlist!

def wave_trans(inimage,wavtranimage):
    iraf.gemini.gmos.gftransform('erg'+inimage, wavtran=wavtranimage)

def qe_corr(inimage,refimage):
    iraf.gemini.gmos.gqecorr(inimage, refimges=refimage, fl_keep=yes)

def extract(inimage

#   start reduction   #
#######################

if setup == "red":

    #load iraf packages
    iraf.load('gemini')
    iraf.load('gmos')

    #start reducing the first flat
    flat1done = False
    while (flat1done == False): 
        run_flats(flat1)
        check1 = os.path.isfile('erg' + flat1base + '_1.fits')
        check2 = os.path.isfile('erg' + flat1base + '_2.fits')
        if (check1 == True) or (check2 == True):
            print('Looks like you did not get all the way through, would you like to start over?')
            print('1=yes')
            print('2=no')
            ans = raw_input('your choice:   ')
            if ans == 1:
                delete(foldername + 'database/' + 'aperg' + flat1base + '_1')
                delete(foldername + 'database/' + 'aperg' + flat1base + '_2')
                delete(foldername + 'erg' + flat1)
                delete(foldername + 'rg' + flat1)
                delete(foldername + 'g' + flat1)
            else: flat1done = True
        else: flat1done = True

    #start reducing the second flat
    flat2done = False
    while (flat2done == False):
        run_flats(flat2)
        check2 = os.path.isfile('erg' + flat2base + '_1.fits')
        if (check1 == True) or (check2 == True):
            print('Looks like you did not get all the way through, would you like to start over?')
            print('1=yes')
            print('2=no')
            ans = raw_input('your choice:   ')
            if ans == 1:
                delete(foldername + 'database/' + 'aperg' + flat2base + '_1')
                delete(foldername + 'database/' + 'aperg' + flat2base + '_2')
                delete(foldername + 'erg' + flat2)
                delete(foldername + 'rg' + flat2)
                delete(foldername + 'g' + flat2)
            else: flat2done = True
        else: flat2done = True

else:
    
    iraf.load('gemini')

    iraf.load('gmos') 

    run_flats(flat1)
    run_flats(flat2)
