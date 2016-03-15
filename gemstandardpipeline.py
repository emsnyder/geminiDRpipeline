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

print('Which star are you reducing?    ')
starname=raw_input('your answer:    ')
print('If this star is an LTT####, now shorten its name to L#### please.')
starnameshort=raw_input('your answer:    ')

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

#~~~~pull out science images and define done parameters

selsci1 = (obsclass == 'partnerCal') & (obj == starname) & (wave == waves[0])
sciimage1 = np.extract(selsci1,files)
n = np.size(sciimage1)
if n <= 1:
    sci1 = sciimage1[0]
    sci1base = sci1.split('.')[0]
    sci1done = 'erg' + sci1base

selsci2 = (obsclass == 'partnerCal') & (obj == starname) & (wave == waves[1])
sciimage2 = np.extract(selsci2,files)
n = np.size(sciimage2)
if n <= 1:
    sci2 = sciimage2[0]
    sci2base = sci2.split('.')[0]
    sci2done = 'erg' + sci2base

#~~~~pull out twilight images and define done parameters

seltwi1 = (obj == 'Twilight') & (wave == waves[0])
twiimage1 = np.extract(seltwi1,files)
n = np.size(twiimage1)
if n <= 1:
    twi1 = twiimage1[0]
    twi1base = twi1.split('.')[0]
    twi1done = 'erg' + twi1base

seltwi2 = (obj == 'Twilight') & (wave == waves[1])
twiimage2 = np.extract(seltwi2,files)
n = np.size(twiimage2)
if n <= 1:
    twi2 = twiimage2[0]
    twi2base = twi2.split('.')[0]
    twi2done = 'erg' + twi2base

#~~~~pull out arcs and define done parameters

selarc1 = (wave == waves[0]) & (obsclass == 'dayCal') & (obj == 'CuAr')
arc1 = np.extract(selarc1, files)
n = np.size(arc1)
if n <= 1:
    arc1 = arc1[0]
    arc1base = arc1.split('.')[0]
    arc1done = 'erg' + arc1base
    arc1wt   = 'terg' + arc1base

selarc2 = (wave == waves[1]) & (obsclass == 'dayCal') & (obj == 'CuAr')
arc2 = np.extract(selarc2, files)
n = np.size(arc2)
if n <= 1:
    arc2 = arc2[0]
    arc2base = arc2.split('.')[0]
    arc2done = 'erg' + arc2base
    arc2wt   = 'terg' + arc2base

#~~~~pull out flats and define done parameters

selflat1 = (wave == waves[0]) & (obsclass == 'partnerCal') & (obj == 'GCALflat')
flat1 = np.extract(selflat1, files)[0]
flat1base = flat1.split('.')[0]
flat1done = 'erg' + flat1base

selflat2 = (wave == waves[1]) & (obsclass == 'partnerCal') & (obj == 'GCALflat')
flat2 = np.extract(selflat2, files)[0]
flat2base = flat2.split('.')[0]
flat2done = 'erg' + flat2base

#~~~~pull out acquisition images

selacq = (obsclass == 'acqCal') & (obj == starname)
acqimages = np.extract(selacq, files)

obsdate = date[0] #date should be same for all files

if obsdate == '2015-11-08': biasim='20151108bias.fits'

#~~~~pull out bias

bias1 = glob.glob('*bias*.fits')[0]
biasbase = bias1.split('.')[0]

#~~~~print out file list
print('--------------------------------------------------')
print('file name            file type    cent. wavelength')
print('--------------------------------------------------')
print(flat1 + '  flat         ' + waves[0])
print(flat2 + '  flat         ' + waves[1])
print(arc1 + '  arc          ' + waves[0])
print(arc2 + '  arc          ' + waves[1])
print(sci1 + '  science      ' + waves[0])
print(sci2 + '  science      ' + waves[1])
print(twi1 + '  twilight     ' + waves[0])
print(twi2 + '  twilight     ' + waves[1])
print(bias1 + '    bias         ---')
print('--------------------------------------------------')
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

#~~~~flat reduction

print('%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flat reduction')
print('%%%%%%%%%%%%%%%%%%%%%%%')
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
                              logfile=starname+'.log', \
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

#~~~~reduce the arcs

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting basic arc reduction')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def reduce_arcs(inimage, refimage, slit):
    iraf.gemini.gmos.gfreduce(inimage, fl_wavtran='no', fl_inter='no', ref=refimage, recenter='no', \
                              trace='no', fl_skysub='no', fl_gscrrej='no', fl_bias='no', fl_over='yes', \
                              fl_fluxcal='no', order='1', weights='none', slits=slit,logfile=starname+'.log')

if os.path.isfile(arc1done+'.fits') == False:
    reduce_arcs(arc1base,flat1done,slit)
    print('>>> arc 1 reduction complete')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)
 
if os.path.isfile(arc2done+'.fits') == False:
    reduce_arcs(arc2base,flat2done,slit)
    print('>>> arc 2 reduction complete')
    time.sleep(2)
else:
    print('>>> arc 2 already done, moving on')
    time.sleep(2)

#~~~~find the wavelength solution for the arcs

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting wavelength solution')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def wave_cal(inimage):
    iraf.gemini.gmos.gswavelength(inimage, fl_inter=yes, nlost=10, coordli='../smalllinelist.dat',logfile=starname+'.log')

if os.path.isfile(arc1wt+'.fits') == False:
    wave_cal('erg' + arc1base)
    print('>>> arc 1 wavelength solution done')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)

if os.path.isfile(arc2wt+'.fits') == False:
    wave_cal('erg' + arc2base)
    print('>>> arc 2 wavelength solution done')
    time.sleep(2)
else:
    print('>>> arc 2 already done, moving on')
    time.sleep(2)

#~~~~wavelength calibration of the arcs

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting wavelength tranformation of arcs')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def wave_trans(inimage,wavtranimage):
    iraf.gemini.gmos.gftransform(inimage, wavtran=wavtranimage,logfile=starname+'.log')

if os.path.isfile(arc1wt+'.fits') == False:
    wave_trans('erg' + arc1base, 'erg' + arc1base)
    print('>>> arc 1 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)

if os.path.isfile(arc2wt+'.fits') == False:
    wave_trans('erg' + arc2base, 'erg' + arc2base)
    print('>>> arc 1 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> arc 1 already done, moving on')
    time.sleep(2)

#~~~~qe correct the flats

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting quantum efficiency correction of flats')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def qe_corr(inimage,refimage):
    iraf.gemini.gmos.gqecorr(inimage, refimages=refimage, fl_keep='yes',logfile=starname+'.log')

if os.path.isfile('qrg' + flat1base + '.fits') == False:
    qe_corr('rg' + flat1base, arc1done)
    print('>>> flat 1 QE correction done')
    time.sleep(2)
else:
    print('>>> flat 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('qrg' + flat2base + '.fits') == False:
    qe_corr('rg' + flat2base, arc2done)
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
    iraf.gemini.gmos.gfextract(inimage,ref=refimage,logfile=starname+'.log')

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

#~~~~reduce the twilights

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting twilight reduction')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def reduce_twis(inimage, refimage, slit):
    iraf.gemini.gmos.gfreduce(inimage, fl_extract='no', fl_wavtran='no', fl_inter='no', ref=refimage, recenter='no', \
                              trace='no', fl_skysub='no', fl_gscrrej='no', fl_bias='yes', fl_over='yes', \
                              order='1', fl_fluxcal='no', weights='none', slits=slit,logfile=starname+'.log',bias=bias1)

if os.path.isfile('rg'+twi1base+'.fits') == False:
    reduce_twis(twi1base,flat1done,slit)
    print('>>> twi 1 reduction complete')
    time.sleep(2)
else:
    print('>>> twi 1 already done, moving on')
    time.sleep(2)
   
if os.path.isfile('rg'+twi2base+'.fits') == False:
    reduce_twis(twi2base,flat2done,slit)
    print('>>> twi 2 reduction complete')
    time.sleep(2)
else:
    print('>>> twi 2 already done, moving on')
    time.sleep(2)

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting quantum efficiency correction of twilights')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def qe_corr2(inimage, refimage, corrimage):
    iraf.gemini.gmos.gqecorr(inimage, refimages=refimage, corrimages=corrimage,logfile=starname+'.log')

if os.path.isfile('qrg' + twi1base + '.fits') == False:
    qe_corr2('rg' + twi1base, arc1done, 'qecorr'+arc1done)
    print('>>> twilight 1 QE correction done')
    time.sleep(2)
else:
    print('>>> twilight 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('qrg' + twi2base + '.fits') == False:
    qe_corr2('rg' + twi2base, arc2done, 'qecorr'+arc2done)
    print('>>> twilight 2 QE correction done')
    time.sleep(2)
else:
    print('>>> twilight 2 already done, moving on')
    time.sleep(2)

#~~~~re-extract the qe-corrected flats

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting extraction of twilights')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def extract(inimage, refimage):
    iraf.gemini.gmos.gfextract(inimage,ref=refimage,logfile=starname+'.log')

if os.path.isfile('eqrg' + twi1base + '.fits') == False:
    extract('qrg' + twi1base, flat1done)
    print('>>> twilight 1 re-extraction done')
    time.sleep(2)
else:
    print('>>> twilight 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('eqrg' + twi2base + '.fits') == False:
    extract('qrg' + twi2base, flat2done)
    print('>>> twilight 2 re-extraction done')
    time.sleep(2)
else:
    print('>>> twilight 2 already done, moving on')
    time.sleep(2)

#~~~~remake the response function

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting response function derivation')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def twi_resp(inimage, outimage, refimage):
    iraf.gemini.gmos.gfresponse(inimage, outimage, sky=refimage, order='95', fl_inter='no', \
                                func='spline3', sample="*")

if os.path.isfile('resp'+twi1base+'.fits') == False:
    twi_resp('eqrg'+twi1base,'resp'+twi1base+'.fits','eqrg'+flat1base)
    print('>>> twi 1 response function complete')
    time.sleep(2)
else:
    print('>>> twi 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('resp'+twi2base+'.fits') == False:
    twi_resp('eqrg'+twi2base,'resp'+twi2base+'.fits','eqrg'+flat2base)
    print('>>> twi 2 response function complete')
    time.sleep(2)
else:
    print('>>> twi 2 already done, moving on')
    time.sleep(2)

response1 = 'resp'+twi1base+'.fits'
response2 = 'resp'+twi2base+'.fits'

#~~~~flatfield, bias/overscan subtract science frames

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flat fielding of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_reduce(inimage,refimage,slit,response):
    iraf.gemini.gmos.gfreduce(inimage, fl_inter='no', fl_addmdf='yes', fl_over='yes', fl_trim='yes', fl_bias='yes', \
                              fl_gscrrej='no', fl_extract='no', fl_wavtran='no', fl_sky='no', fl_fluxcal='no', \
                              slits=slit, trace='no', verb='yes', refer=refimage, response=response, weights='none', \
                              fl_qecorr='no', fl_crspec='no',logfile=starname+'.log', bias=bias1)

if os.path.isfile('rg' + sci1base + '.fits') == False:
    sci_reduce(sci1base,'eqrg'+flat1base,slit,response1)
    print('>>> sci 1 flat fielding done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('rg' + sci2base + '.fits') == False:
    sci_reduce(sci2base,'eqrg'+flat2base,slit,response2)
    print('>>> sci 2 flat fielding done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

#~~~~ remove cosmic rays in science frames

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting CR rejection of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_cr(inimage,outimage):
    iraf.gemini.gemcrspec(inimage,outimage,sigclip='3',fl_vardq='yes')

if os.path.isfile('xrg' + sci1base + '.fits') == False:
    sci_cr('rg'+sci1base,'xrg'+sci1base)
    print('>>> sci 1 CR rejection done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('xrg' + sci2base + '.fits') == False:
    sci_cr('rg'+sci2base,'xrg'+sci2base)
    print('>>> sci 2 CR rejection done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

#~~~~QE correct science frames

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting quantum efficiency correction of science frames')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def qe_corr2(inimage, refimage, corrimage):
    iraf.gemini.gmos.gqecorr(inimage, refimages=refimage, corrimages=corrimage,logfile=starname+'.log')

if os.path.isfile('qxrg' + sci1base + '.fits') == False:
    qe_corr2('xrg' + sci1base, arc1done, 'qecorr'+arc1done)
    print('>>> science 1 QE correction done')
    time.sleep(2)
else:
    print('>>> science 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('qxrg' + sci2base + '.fits') == False:
    qe_corr2('xrg' + sci2base, arc2done, 'qecorr'+arc2done)
    print('>>> science 2 QE correction done')
    time.sleep(2)
else:
    print('>>> science 2 already done, moving on')
    time.sleep(2)

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting extraction of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def extract(inimage, refimage):
    iraf.gemini.gmos.gfextract(inimage,ref=refimage,logfile=starname+'.log')

if os.path.isfile('eqxrg' + sci1base + '.fits') == False:
    extract('qxrg' + sci1base, flat1done)
    print('>>> science 1 extraction done')
    time.sleep(2)
else:
    print('>>> science 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('eqxrg' + sci2base + '.fits') == False:
    extract('qxrg' + sci2base, flat2done)
    print('>>> science 2 extraction done')
    time.sleep(2)
else:
    print('>>> science 2 already done, moving on')
    time.sleep(2)

#~~~~wavelength calibrate science frames!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting wavelength calibration of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_wave(inimage,refimage):
    iraf.gemini.gmos.gftransform(inimage, wavtran=refimage,logfile=starname+'.log')

if os.path.isfile('teqxrg' + sci1base + '.fits') == False:
    sci_wave('eqxrg' + sci1base,arc1done)
    print('>>> sci 1 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('teqxrg' + sci2base + '.fits') == False:
    sci_wave('eqxrg' + sci2base,arc2done)
    print('>>> sci 2 wavelength transformation done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

quit

#~~~~subtract the sky!

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting sky subtraction of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_sky(inimage):
    iraf.gemini.gmos.gfskysub(inimage, fl_inter='yes', verb='yes', weight='none',sepslits='yes',logfile=starname+'.log')

if os.path.isfile('steqxrg' + sci1base + '.fits') == False:
    sci_sky('teqxrg' + sci1base)
    print('>>> sci 1 sky subtraction done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('steqxrg' + sci2base + '.fits') == False:
    sci_sky('teqxrg' + sci2base)
    print('>>> sci 2 sky subtraction  done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

#~~~~Make the Sensitivity Curve
# 1--sum all spectra

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting sensitivity curve creation')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sum_aps(inimage, outimage):
    iraf.gemini.gmos.gfapsum(inimage, outimages=outimage, combine="sum", reject="pclip", fl_inter='no',logfile=starname+'.log')

if os.path.isfile('apsum_steqxrg' + sci1base + '.fits') == False:
    sum_aps('steqxrg' + sci1base, 'apsum_steqxrg' + sci1base + '.fits')
    print('>>> sci 1 aperture summation done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('apsum_steqxrg' + sci2base + '.fits') == False:
    sum_aps('steqxrg' + sci2base, 'apsum_steqxrg' + sci2base + '.fits')
    print('>>> sci 2 aperture summation done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

# 2--make the curve

def sci_standard(inimage, sfile, sfunction, star):
    iraf.gemini.gmos.gsstandard(inimage, sfile=sfile, sfunctio=sfunction, starname=star, fl_inte='no', \
                                observa="Gemini-South", functio="chebyshev", order="4", caldir="onedstds$ctionewcal/",logfile=starname+'.log')

if os.path.isfile('sfunction_steqxrg' + sci1base + '.fits') == False:
    sci_standard('apsum_steqxrg' + sci1base + '.fits','sfile_steqxrg' + sci1base + '.fits','sfunction_steqxrg' + sci1base + '.fits',starnameshort)
    print('>>> sci 1 sensitivity curve created')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('sfunction_steqxrg' + sci2base + '.fits') == False:
    sci_standard('apsum_steqxrg' + sci2base + '.fits','sfile_steqxrg' + sci2base + '.fits','sfunction_steqxrg' + sci2base + '.fits',starnameshort)
    print('>>> sci 2 sensitivity curve created')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

#~~~~Flux calibrate the spectra

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting flux calibration of science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_calibrate(inimage, sfunc):
    iraf.gemini.gmos.gscalibrate(inimage,sfuncti=sfunc,observa="Gemini-South",fluxscal='1',logfile=starname+'.log')

if os.path.isfile('capsum_steqxrg' + sci1base + '.fits') == False:
    sci_calibrate('apsum_steqxrg' + sci1base + '.fits','sfunction_steqxrg' + sci1base + '.fits')
    print('>>> sci 1 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('capsum_steqxrg' + sci2base + '.fits') == False:
    sci_calibrate('apsum_steqxrg' + sci2base + '.fits','sfunction_steqxrg' + sci2base + '.fits')
    print('>>> sci 2 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)

#~~~~create data cubes

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('starting data cube creation from science')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
time.sleep(2)

def sci_cube(inimage,sam):
    iraf.gemini.gmos.gfcube(inimage, ssample=sam,logfile=starname+'.log')

if os.path.isfile('dsteqxrg' + sci11base + '.fits') == False:
    sci_cube('csteqxrg' + sci11base,'0.2')
    print('>>> sci 1 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 1 already done, moving on')
    time.sleep(2)

if os.path.isfile('dcsteqxrg' + sci21base + '.fits') == False:
    sci_cube('csteqxrg' + sci21base,'0.2')
    print('>>> sci 2 flux calibration done')
    time.sleep(2)
else:
    print('>>> sci 2 already done, moving on')
    time.sleep(2)
