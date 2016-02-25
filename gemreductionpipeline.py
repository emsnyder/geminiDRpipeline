import numpy as np
import pylab
import pyfits
import sys
import glob
import os
from pyraf import iraf
import glob
import time 
from astropy.io import fits as fits

print '----------------------------------------------'
print 'Welcome to the Gemini data reduction pipeline!'
print '----------------------------------------------'
print ' '
print ' '

############################
#find where user is working#
############################

mainpath = os.getcwd()

foldername=raw_input('What is the folder you are working in?    ')

if mainpath+'/' != foldername:
    print('you are not in the correct directory! quitting')
    sys.exit()

##############################################
#sort all fits files in the working directory#
##############################################

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

#~~~~select either the red or blue setups

if setup == 'blue':
    waves = np.array(['480.0', '490.0'])
    slit = 'red'
    vers = '1'
if setup == 'red':
    waves = np.array(['635.0', '640.0'])
    slit = 'both'
    vers = '*'

#~~~~pull out acquisition images

selacq = (obsclass == 'acq')
acqimages = np.extract(selacq, files)

#~~~~pull out science images and define done parameters

selsci1 = (obsclass == 'science') & (wave == waves[0])
sciimage1 = np.extract(selsci1,files)
n = np.size(sciimage1)
if n <= 1:
    sci11 = sciimage1[0]
    sci11base = sci11.split('.')[0]
    sci11done = 'erg' + sci11
else:
    sci11 = sciimage1[0]
    sci12 = sciimage1[1]
    sci11base = sci11.split('.')[0]
    sci12base = sci12.split('.')[0]
    sci11done = 'erg' + sci11
    sci12done = 'erg' + sci12

selsci2 = (obsclass == 'science') & (wave == waves[1])
sciimage2 = np.extract(selsci2,files)
if n <= 1:
    sci21 = sciimage2[0]
    sci21base = sci21.split('.')[0]
    sci21done = 'erg' + sci21
else:
    sci21 = sciimage2[0]
    sci22 = sciimage2[1]
    sci21base = sci21.split('.')[0]
    sci22base = sci22.split('.')[0]
    sci21done = 'erg' + sci21
    sci22done = 'erg' + sci22

#~~~~pull out arcs and define done parameters

selarc1 = (wave == waves[0]) & (obstype == 'ARC')
arc1 = np.extract(selarc1, files)
n = np.size(arc1)
if n <= 1:
    arc11 = arc1[0]
    arc11base = arc11.split('.')[0]
    arc11done = 'erg' + arc11
    arc11wt   = 'terg' + arc11
else:
    arc11 = arc1[0]
    arc12 = arc1[1]
    arc11base = arc11.split('.')[0]
    arc12base = arc12.split('.')[0]
    arc11done = 'erg' + arc11
    arc12done = 'erg' + arc12
    arc11wt   = 'terg' + arc11 #see if wavetran applied to arc
    arc12wt   = 'terg' + arc12

selarc2 = (wave == waves[1]) & (obstype == 'ARC')
arc2 = np.extract(selarc2, files)
n = np.size(arc2)
if n <= 1:
    arc21 = arc2[0]
    arc21base = arc21.split('.')[0]
    arc21done = 'erg' + arc21
    arc21wt   = 'terg' + arc21
else:
    arc21 = arc2[0]
    arc22 = arc2[1]
    arc21base = arc21.split('.')[0]
    arc22base = arc22.split('.')[0]
    arc21done = 'erg' + arc21
    arc22done = 'erg' + arc22
    arc21wt   = 'terg' + arc21
    arc22wt   = 'terg' + arc22

#~~~~pull out flats and define done parameters

selflat1 = (wave == waves[0]) & (obstype == 'FLAT')
flat1 = np.extract(selflat1, files)[0]
flat1base = flat1.split('.')[0]
flat1done = 'erg' + flat1
selflat2 = (wave == waves[1]) & (obstype == 'FLAT')
flat2 = np.extract(selflat2, files)[0]
flat2base = flat2.split('.')[0]
flat2done = 'erg' + flat2

#~~~~pull out date

obsdate = date[0] #date should be same for all files

if obsdate == '2015-11-08': biasim='20151108bias.fits'

#~~~~pull out bias

bias1 = glob.glob('*bias*.fits')[0]
biasbase = bias1.split('.')[0]

#~~~~pull out standard star info

#if date = 2015 and setup = blue then standard star folder = ??
 
#~~~~pull out galaxy name

selname = (obsclass == 'science') & (wave == waves[0])
galname = np.extract(selname, obj)

#~~~~print out file list
print('file name            file type    cent. wavelength')
print(flat1 + '  flat         ' + waves[0])
print(flat2 + '  flat         ' + waves[1])
print(arc11 + '  arc          ' + waves[0])
try:
    arc12
except NameError:
    pass 
else:
    print(arc12 + '  arc       ' + waves[0])
print(arc21 + '  arc          ' + waves[1])
try:
    arc12
except NameError:
    pass 
else:
    print(arc22 + ' arc         ' + waves[1])

print(sci11 + '  science      ' + waves[0])
try:
    sci12
except NameError:
    pass 
else:
    print(sci12 + '  science      ' + waves[0])
print(sci21 + '  science      ' + waves[1])
try:
    sci22
except NameError:
    pass
else:
    print(sci22 + '  science      ' + waves[1])
print(bias1 + '    bias         ---')

print('')
print('')
print('Do these look correct to you?')
print('1=yes')
print('0=no')
ans = raw_input('your choice:   ')

if ans == '0':
    print('Whoa, a bug already?! Shoot Elaine an email.')
    sys.exit()
else:
    print('Alright, on to the reduction...')
    time.sleep(3)

###############################
#define all the IRAF functions#
###############################

def run_gfreduce(inimage,biasimage,refimage,slit,cr,wt,ss,inter,os,tr,bi,ex,fl,ap,weight,trac,rec,orde):
    iraf.gemini.gmos.gfreduce(inimage, \
                              fl_inter=inter, \
                              fl_over=os, \
                              fl_trim=tr, \
                              fl_bias=bi, \
                              fl_flux=fl, \
                              fl_extract=ex, \
                              fl_gscrrej=cr, \
                              fl_wavtran=wt,
                              fl_gsappwave=ap, \
                              fl_skysub=ss, \
                              slits=slit, \
                              rawpath=foldername, \
                              weights=weight, \
                              trace=trac, \
                              recenter=rec, \
                              order=orde, \
                              logfile=galname[0]+'.log', \
                              bias=biasimage)

def plot_spaxels(inimage,verinfo):
    iraf.gemini.gmos.gfdisplay(inimage, ver=verinfo)

def wave_cal(inimage):
    iraf.gemini.gmos.gswavelength(inimage, fl_inter=yes, nlost=10)
#add coordlist!

def wave_trans(inimage,wavtranimage):
    iraf.gemini.gmos.gftransform('erg'+inimage, wavtran=wavtranimage)

def qe_corr(inimage,refimage):
    iraf.gemini.gmos.gqecorr(inimage, refimges=refimage, fl_keep=yes)


#################
#start reduction#
#################

#~~~~load iraf packages
iraf.load('gemini')
time.sleep(3)
iraf.load('gmos')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flat reduction and fiber identification')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(3)



#~~~~flat 1 reduction

#if os.path.isfile(flat1done) == False: #if the output file isn't done yet, do this step

#    flat1done = False
#    while (flat1done == False): 
#        run_gfreduce(flat1,bias1,'',slit,'no','no','no','yes','yes','yes','yes','yes','no','no','none','no','no','INDEF')
#        check1 = os.path.isfile('erg' + flat1base + '_1.fits')
#        check2 = os.path.isfile('erg' + flat1base + '_2.fits')
#        if (check1 == True) or (check2 == True):
#            print('Looks like you did not get all the way through, would you like to start over?')
#            print('1=yes')
#            print('0=no')
#            ans = raw_input('your choice:   ')
#            if ans == 1:
#                delete(foldername + 'database/' + 'aperg' + flat1base + '_1')
#                delete(foldername + 'database/' + 'aperg' + flat1base + '_2')
#                delete(foldername + 'erg' + flat1base + '_1.fits')
#                delete(foldername + 'erg' + flat1base + '_2.fits')
#                delete(foldername + 'erg' + flat1)
#                delete(foldername + 'rg' + flat1)
#                delete(foldername + 'g' + flat1)
#            else: flat1done = True
#        else: flat1done = True

#else:
#    print('the first flat is already done, moving on')

#~~~~flat 2 reduction

#if os.path.isfile(flat2done) == False:

#    flat2done = False
#    while (flat2done == False):
#        run_gfreduce(flat2,bias1,'',slit,'no','no','no','yes','yes','yes','yes','yes','no','no','none','no','no','')
#        check1 = os.path.isfile('erg' + flat2base + '_1.fits')
#        check2 = os.path.isfile('erg' + flat2base + '_2.fits')
#        if (check1 == True) or (check2 == True):
#            print('Looks like you did not get all the way through, would you like to start over?')
#            print('1=yes')
#            print('0=no')
#            ans = raw_input('your choice:   ')
#            if ans == 1:
#                delete(foldername + 'database/' + 'aperg' + flat2base + '_1')
#                delete(foldername + 'database/' + 'aperg' + flat2base + '_2')
#                delete(foldername + 'erg' + flat2base + '_1.fits')
#                delete(foldername + 'erg' + flat2base + '_2.fits')
#                delete(foldername + 'erg' + flat2)
#                delete(foldername + 'rg' + flat2)
#                delete(foldername + 'g' + flat2)
#            else: flat2done = True
#        else: flat2done = True

#else:
#    print('the second flat is already done, moving on')


if os.path.isfile(flat1done) == False: #if the output file isn't done yet, do this step
    run_gfreduce(flat1,bias1,'',slit,'no','no','no','yes','yes','yes','yes','yes','no','no','none','yes','yes','default')

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('flat 1 complete, on to flat 2')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(3)

if os.path.isfile(flat2done) == False:
    run_gfreduce(flat2,bias1,'',slit,'no','no','no','yes','yes','yes','yes','yes','no','no','none','yes','yes','default')

print('%%%%%%%%%%%%%%%')
print('flat 2 complete')
print('%%%%%%%%%%%%%%%')
time.sleep(2)
print('%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting arc 1 reduction')
print('%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

#~~~~arc reduction

if os.path.isfile(arc11done) == False:
    run_gfreduce(arc11,bias1,flat1done,slit,'no','no','no','no','yes','yes','yes','yes','no','no','none','no','no','1')

try:
    arc12
except NameError:
    pass
else:
    if os.path.isfile(arc12done) == False:
        run_gfreduce(arc12,bias1,flat1done,slit,'no','no','no','no','yes','no','no','yes','no','no','none','no','no','1')

if os.path.isfile(arc21done) == False:
    run_gfreduce(arc21,bias1,flat2done,slit,'no','no','no','no','yes','no','no','yes','no','no','none','no','no','1')

try:
    arc22
except NameError:
    pass
else:
    if os.path.isfile(arc22done) == False:
        run_gfreduce(arc22,bias1,flat2done,slit,'no','no','no','no','yes','no','no','yes','no','no','none','no','no','1')

#~~~~wavelength solution

