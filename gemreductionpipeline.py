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

print "    .         *             .     . +    .     .          .       "
print "             .      .    .             .         \|/              "
print "        *                       +          `    - * -         `   "
print "  `         .     `       +   .    *          .  /|\              "
print "                   +                 .     .             +      . "
print "  .      .  *           *     +     .     .        *              "
print "    .      .     .          .   .           .        .        *   "
print "          ----------------------------------------------          "
print "  *       Welcome to the Gemini data reduction pipeline!   `      "
print "          ----------------------------------------------          "
print "    .    _     *       \|/   .       .      -*-              +    "
print "      .` \\`.     +     -*-     *   .         `       .   *       "
print "   .  |__''_|  .       /|\ +         .    +       .           |   "
print "      |     | .                                        .     -o-  "
print "      |     |           `  .        `           ,'`.   *   .  |   "
print "    _.'-----'-._     *         +    ,'`. ,'`. ,'    `.            "
print "  /             \__.__.--._________'    `.   `.       `._________ "



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
    sci11done = 'erg' + sci11base
else:
    sci11 = sciimage1[0]
    sci12 = sciimage1[1]
    sci11base = sci11.split('.')[0]
    sci12base = sci12.split('.')[0]
    sci11done = 'erg' + sci11base
    sci12done = 'erg' + sci12base

selsci2 = (obsclass == 'science') & (wave == waves[1])
sciimage2 = np.extract(selsci2,files)
n = np.size(sciimage2)
if n <= 1:
    sci21 = sciimage2[0]
    sci21base = sci21.split('.')[0]
    sci21done = 'erg' + sci21base
else:
    sci21 = sciimage2[0]
    sci22 = sciimage2[1]
    sci21base = sci21.split('.')[0]
    sci22base = sci22.split('.')[0]
    sci21done = 'erg' + sci21base
    sci22done = 'erg' + sci22base

#~~~~pull out arcs and define done parameters

selarc1 = (wave == waves[0]) & (obstype == 'ARC')
arc1 = np.extract(selarc1, files)
n = np.size(arc1)
if n <= 1:
    arc11 = arc1[0]
    arc11base = arc11.split('.')[0]
    arc11done = 'erg' + arc11base
    arc11wt   = 'terg' + arc11base
else:
    arc11 = arc1[0]
    arc12 = arc1[1]
    arc11base = arc11.split('.')[0]
    arc12base = arc12.split('.')[0]
    arc11done = 'erg' + arc11base
    arc12done = 'erg' + arc12base
    arc11wt   = 'terg' + arc11base #see if wavetran applied to arc
    arc12wt   = 'terg' + arc12base

selarc2 = (wave == waves[1]) & (obstype == 'ARC')
arc2 = np.extract(selarc2, files)
n = np.size(arc2)
if n <= 1:
    arc21 = arc2[0]
    arc21base = arc21.split('.')[0]
    arc21done = 'erg' + arc21base
    arc21wt   = 'terg' + arc21base
else:
    arc21 = arc2[0]
    arc22 = arc2[1]
    arc21base = arc21.split('.')[0]
    arc22base = arc22.split('.')[0]
    arc21done = 'erg' + arc21base
    arc22done = 'erg' + arc22base
    arc21wt   = 'terg' + arc21
    arc22wt   = 'terg' + arc22

#~~~~pull out flats and define done parameters

selflat1 = (wave == waves[0]) & (obstype == 'FLAT')
flat1 = np.extract(selflat1, files)[0]
flat1base = flat1.split('.')[0]
flat1done = 'erg' + flat1base
selflat2 = (wave == waves[1]) & (obstype == 'FLAT')
flat2 = np.extract(selflat2, files)[0]
flat2base = flat2.split('.')[0]
flat2done = 'erg' + flat2base

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
    arc22
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
    time.sleep(2)

#################
#start reduction#
#################

#~~~~load iraf packages
iraf.load('gemini')
time.sleep(2)
iraf.load('gmos')

#~~~~flat 1 reduction

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flat reduction and fiber identification')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def run_gfreduce(inimage,biasimage,refimage,slit,cr,wt,ss,inter,os,tr,bi,ex,fl,ap,weight,trac,rec,orde):
    iraf.gemini.gmos.gfreduce(inimage, \
                              fl_inter=inter, \
                              fl_over=os, \
                              fl_trim=tr, \
                              fl_bias=bi, \
                              fl_fluxcal=fl, \
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

if os.path.isfile(flat1done+'.fits') == False: #if the output file isn't done yet, do this step
    run_gfreduce(flat1,bias1,'',slit,'no','no','no','yes','yes','yes','yes','yes','no','no','none','yes','yes','default')
    print('>>> flat 1 complete, on to flat 2')
    time.sleep(2)
else:
    print('>>> flat 1 already done, moving on')
    time.sleep(2)

if os.path.isfile(flat2done+'.fits') == False:
    run_gfreduce(flat2,bias1,'',slit,'no','no','no','yes','yes','yes','yes','yes','no','no','none','yes','yes','default')
    print('>>> flat 2 complete')
    time.sleep(2)
else:
    print('>>> flat 2 already done, moving on')
    time.sleep(2)

#~~~~basic arc reduction

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting basic arc reduction')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def reduce_arcs(inimage, refimage, slit):
    iraf.gemini.gmos.gfreduce(inimage, fl_wavtran='no', fl_inter='no', ref=refimage, recenter='no', \
                              trace='no', fl_skysub='no', fl_gscrrej='no', fl_bias='no', fl_over='yes', \
                              order='1', weights='none', slits=slit)

if os.path.isfile(arc11done+'.fits') == False:
    reduce_arcs(arc11base,flat1done,slit)
    print('>>> arc 1 reduction complete')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)
   
try:
    arc12
except NameError:
    pass
else:
    if os.path.isfile(arc12done+'.fits') == False:
        reduce_arcs(arc12base,flat1done,slit)
        print('>>> arc 12 reduction complete')
        time.sleep(2)
    else:
        print('>>> arc 12 already done, moving on')
        time.sleep(2)

if os.path.isfile(arc21done+'.fits') == False:
    reduce_arcs(arc21base,flat2done,slit)
    print('>>> arc 2 reduction complete')
    time.sleep(2)
else:
    print('>>> arc 2 already done, moving on')
    time.sleep(2)

try:
    arc22
except NameError:
    pass
else:
    if os.path.isfile(arc22done+'.fits') == False:
        reduce_arcs(arc22base,flat2done,slit)
        print('>>> arc 22 reduction complete')
        time.sleep(2)
    else:
        print('>>> arc 22 already done, moving on')
        time.sleep(2)

#~~~~find the wavelength solution for the arcs

coordlistloc = raw_input('Where is your line list?     ')

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting wavelength solution')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def wave_cal(inimage):
    iraf.gemini.gmos.gswavelength(inimage, fl_inter=yes, nlost=10,coordli=coordlistloc)

if os.path.isfile(arc11wt+'.fits') == False:
    wave_cal('erg' + arc11base)
    print('>>> arc 1 wavelength solution done')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)

try:
    arc12
except NameError:
    pass
else:
    if os.path.isfile(arc12wt+'.fits') == False:
        wave_cal('erg' + arc12base)
        print('>>> arc 12 wavelength solution done')
        time.sleep(2)
    else:
        print('>>> arc 12 already done, moving on')
        time.sleep(2)

if os.path.isfile(arc21wt+'.fits') == False:
    wave_cal('erg' + arc21base)
    print('>>> arc 2 wavelength solution done')
    time.sleep(2)
else:
    print('>>> arc 2 already done, moving on')
    time.sleep(2)

try:
    arc22
except NameError:
    pass
else:
    if os.path.isfile(arc22wt+'.fits') == False:
        wave_cal('erg' + arc22base)
        print('>>> arc 22 wavelength solution done')
        time.sleep(2)
    else:
        print('>>> arc 22 already done, moving on')
        time.sleep(2)

#~~~~wavelength calibration of the arcs

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting wavelength tranformation of arcs')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def wave_trans(inimage,wavtranimage):
    iraf.gemini.gmos.gftransform(inimage, wavtran=wavtranimage)

if os.path.isfile(arc11wt+'.fits') == False:
    wave_trans('erg' + arc11base, 'erg' + arc11base)
    print('>>> arc 1 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)

try:
    arc12
except NameError:
    pass
else:
    if os.path.isfile(arc12wt+'.fits') == False:
        wave_trans('erg' + arc12base,'erg' + arc12base)
        print('>>> arc 12 wavelength transformation done')
        time.sleep(2)
    else:
        print('>>> arc 12 already done, moving on')
        time.sleep(2)


if os.path.isfile(arc21wt+'.fits') == False:
    wave_trans('erg' + arc21base,'erg' + arc21base)
    print('>>> arc 2 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> arc 2 already done, moving on')
    time.sleep(2)

try:
    arc22
except NameError:
    pass
else:
    if os.path.isfile(arc22wt+'.fits') == False:
        wave_trans('erg' + arc22base,'erg' + arc22base)
        print('>>> arc 22 wavelength transformation done')
        time.sleep(2)
    else:
        print('>>> arc 22 already done, moving on')
        time.sleep(2)

#~~~~qe correct the flats

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting quantum efficiency correction of flats')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def qe_corr(inimage,refimage):
    iraf.gemini.gmos.gqecorr(inimage, refimages=refimage, fl_keep=yes)

if os.path.isfile('qrg' + flat1base + '.fits') == False:
    qe_corr('rg' + flat1base, arc11done)
    print('>>> flat 1 QE correction done')
    time.sleep(2)
else:
    print('>>> flat 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('qrg' + flat2base + '.fits') == False:
    qe_corr('rg' + flat2base, arc21done)
    print('>>> flat 2 QE correction done')
    time.sleep(2)
else:
    print('>>> flat 2 already done, moving on')
    time.sleep(2)

#~~~~re-extract the qe-corrected flats

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting re-extraction of flats')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def extract(inimage, refimage):
    iraf.gemini.gmos.gfextract(inimage,ref=refimage)

if os.path.isfile('eqrg' + flat1base + '.fits') == False:
    extract('qrg' + flat1base, flat1done)
    print('>>> flat 1 re-extraction done')
    time.sleep(2)
else:
    print('>>> flat 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('eqrg' + flat2base + '.fits') == False:
    extract('qrg' + flat2base, flat2done)
    print('>>> flat 2 re-extraction done')
    time.sleep(2)
else:
    print('>>> flat 2 already done, moving on')
    time.sleep(2)

#~~~~bias and overscan subtract the science data

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting bias and overscan subtraction of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_red_bias(inimage,slit):
    iraf.gemini.gmos.gfreduce(inimage, slits=slit, fl_inter='no', fl_over='yes', fl_trim='yes', fl_bias='yes', \
                              fl_flux='no', fl_gscrrej='no', fl_extract='no', fl_gsappwave='no', fl_wavtran='no', \
                              fl_skysub='no', weights='none', bias=bias1)

if os.path.isfile('rg' + sci11base + '.fits') == False:
    sci_red_bias(sci11base,slit)
    print('>>> sci 1 basic reduction done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('rg' + sci12base + '.fits') == False:
        sci_red_bias(sci12base,slit)
        print('>>> sci 11 basic reduction done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('rg' + sci21base + '.fits') == False:
    sci_red_bias(sci21base,slit)
    print('>>> sci 2 basic reduction done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('rg' + sci22base + '.fits') == False:
        sci_red_bias(sci22base,slit)
        print('>>> sci 22 basic reduction done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~remove scattered light

#~~~~remove the cosmic rays

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting cosmic ray rejection on science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def gem_crspec(inimage, outimage):
    iraf.gemini.gmos.gemcrspec(inimage, outiamge)

if os.path.isfile('xrg' + sci11base + '.fits') == False:
    gem_crspec('rg' + sci11base,'xrg' + sci22base)
    print('>>> sci 1 cosmic ray rejection done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('xrg' + sci12base + '.fits') == False:
        gem_crspec('rg' + sci12base,'xrg' + sci22base)
        print('>>> sci 11 cosmic ray rejection done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('xrg' + sci21base + '.fits') == False:
    gem_crspec('rg' + sci21base,'xrg' + sci22base)
    print('>>> sci 2 cosmic ray rejection done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('xrg' + sci22base + '.fits') == False:
        gem_crspec('rg' + sci22base,'xrg' + sci22base)
        print('>>> sci 22 cosmic ray rejection done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~qe correct the science data

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting QE corretion on science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sciqe_corr(inimage,refimage,corr):
    iraf.gemini.gmos.gqecorr(inimage, refimages=refimage, corrimages=corr)

if os.path.isfile('qxrg' + sci11base + '.fits') == False:
    sciqe_corr('xrg' + sci11base,arc11done,'qecorr'+arc11done+'.fits')
    print('>>> sci 1 QE correction done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('qxrg' + sci12base + '.fits') == False:
        sciqe_corr('xrg' + sci12base,arc11done,'qecorr'+arc12done+'.fits')
        print('>>> sci 11 QE correction done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('qxrg' + sci21base + '.fits') == False:
    sciqe_corr('xrg' + sci21base,arc21done,'qecorr'+arc21done+'.fits')
    print('>>> sci 2 QE correction  done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('qxrg' + sci22base + '.fits') == False:
        sciqe_corr('xrg' + sci22base,arc21done,'qecorr'+arc21done+'.fits')
        print('>>> sci 22 QE correction  done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~Flat field and extract spectra!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flat fielding and extraction of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_extract(inimage,refimage,slit,response):
    iraf.gemini.gmos.gfreduce(inimage, fl_inter='no', fl_addmdf='no', fl_over='no', fl_trim='no', fl_bias='no', \
                              fl_gscrrej='no', fl_extract='yes', fl_wavtran='no', fl_sky='no', fl_flux='no', \
                              slits=slit, trace='no', verb='yes', refer=refimage, response=response, weights='none')

if os.path.isfile('eqxrg' + sci11base + '.fits') == False:
    sci_extract('qxrg' + sci11base,'eqrg'+flat1base,slit,response1)
    print('>>> sci 1 extraction done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('eqxrg' + sci12base + '.fits') == False:
        sci_extract('qxrg' + sci12base,'eqrg'+flat1base,slit,response1)
        print('>>> sci 11 extraction done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('eqxrg' + sci21base + '.fits') == False:
    sci_extract('qxrg' + sci21base,'eqrg'+flat2base,slit,response2)
    print('>>> sci 2 extraction  done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('eqxrg' + sci22base + '.fits') == False:
        sci_extract('qxrg' + sci22base,'eqrg'+flat2base,slit,response2)
        print('>>> sci 22 extraction  done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~wavelength calibrate science frames!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting wavelength calibration of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_wave(inimage,refimage):
    iraf.gemini.gmos.gftransform(inimage, wavtran=refimage)

if os.path.isfile('teqxrg' + sci11base + '.fits') == False:
    sci_wave('eqxrg' + sci11base,arc11done)
    print('>>> sci 1 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('teqxrg' + sci12base + '.fits') == False:
        sci_wave('eqxrg' + sci12base,arc11done)
        print('>>> sci 11 wavelength transformation done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('teqxrg' + sci21base + '.fits') == False:
    sci_wave('eqxrg' + sci21base,arc21done)
    print('>>> sci 2 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('teqxrg' + sci22base + '.fits') == False:
        sci_wave('eqxrg' + sci22base,arc21done)
        print('>>> sci 22 wavelength transformation done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~subtract the sky!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting sky subtraction of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_sky(inimage,slit):
    iraf.gemini.gmos.gfreduce(inimage, fl_inter='no', fl_addmdf='no', fl_over='no', fl_trim='no', fl_bias='no', \
                              fl_gscrrej='no', fl_extract='no', fl_wavtran='no', fl_sky='yes', fl_flux='no', \
                              slits=slit, verb='yes', weights='none')

if os.path.isfile('steqxrg' + sci11base + '.fits') == False:
    sci_sky('teqxrg' + sci11base)
    print('>>> sci 1 sky subtraction done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('steqxrg' + sci12base + '.fits') == False:
        sci_sky('teqxrg' + sci12base)
        print('>>> sci 11 sky subtraction done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('steqxrg' + sci21base + '.fits') == False:
    sci_sky('teqxrg' + sci21base)
    print('>>> sci 2 sky subtraction  done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('steqxrg' + sci22base + '.fits') == False:
        sci_sky('teqxrg' + sci22base)
        print('>>> sci 22 sky subtraction  done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~flux calibration!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flux calibration of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_flux(inimage,sfunc):
    iraf.gemini.gmos.gscalibrate(inimage, sfunctio=sfunc, observa="Gemini-South", fluxscal=1)

if os.path.isfile('csteqxrg' + sci11base + '.fits') == False:
    sci_flux('steqxrg' + sci11base,sfunc)
    print('>>> sci 1 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('csteqxrg' + sci12base + '.fits') == False:
        sci_flux('steqxrg' + sci12base,sfunc)
        print('>>> sci 11 flux calibration done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('csteqxrg' + sci21base + '.fits') == False:
    sci_flux('steqxrg' + sci21base,sfunc)
    print('>>> sci 2 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('csteqxrg' + sci22base + '.fits') == False:
        sci_fkux('steqxrg' + sci22base,sfunc)
        print('>>> sci 22 flux calibration  done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

#~~~~create data cubes!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting data cube creation from science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_cube(inimage,sam):
    iraf.gemini.gmos.gfcube(inimage, ssample=sam)

if os.path.isfile('dcsteqxrg' + sci11base + '.fits') == False:
    sci_cube('csteqxrg' + sci11base,sam)
    print('>>> sci 1 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

try:
    sci12
except NameError:
    pass
else:
    if os.path.isfile('dcsteqxrg' + sci12base + '.fits') == False:
        sci_cube('csteqxrg' + sci12base,sam)
        print('>>> sci 11 flux calibration done')
        time.sleep(2)
    else:
        print('>>> sci 11 already done, moving on')
        time.sleep(2)

if os.path.isfile('dcsteqxrg' + sci21base + '.fits') == False:
    sci_cube('csteqxrg' + sci21base,sam)
    print('>>> sci 2 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

try:
    sci22
except NameError:
    pass
else:
    if os.path.isfile('dcsteqxrg' + sci22base + '.fits') == False:
        sci_cube('csteqxrg' + sci22base,sam)
        print('>>> sci 22 flux calibration  done')
        time.sleep(2)
    else:
        print('>>> sci 22 already done, moving on')
        time.sleep(2)

