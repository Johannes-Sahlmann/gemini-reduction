"""
# J. Sahlmann, STScI/AURA, 2016-12-12 -- 2017-02-07
# https://archive.gemini.edu/help/api.html

"""

from __future__ import print_function

# Load the required packages
import os, sys
import numpy as np
from astropy import units as u
from astropy.table import Table, Column

import aplpy
import pylab as pl

import urllib
import json
# from astropy.io import fits
from astropy.time import Time, TimeDelta


from pyraf import iraf
from iraf import gemini
from pyraf.iraf import gmos

iraf.unlearn(gemini)
iraf.unlearn(gmos)

################################################################################

def make_dir(myDir):	
    if not os.path.exists(myDir):
        os.makedirs(myDir)

def get_time_window(tdate,timeDeltaDays):
    deltaDate = TimeDelta(timeDeltaDays)
    t1 = (tdate - deltaDate).iso.split(' ')[0].replace('-','')
    t2 = (tdate + deltaDate).iso.split(' ')[0].replace('-','')
    return t1,t2

def get_json_summary(url):
    uu = urllib.urlopen(url);    
    jsondoc = uu.read();    
    uu.close();    
    files = json.loads(jsondoc);
    return files
    
################################################################################

class JsonSummary(object):
    """a structure class for jason summary from gemini archive"""
    def __init__(self, url=None):
        self.url = url
        self.files = get_json_summary(url)
        print('Gemini archive contains %d files with this query: %s'% (len(self.files),url))

    def print_full_info(self):        
        for f in self.files:
            for key in f.keys():
                print(f[key],end='  ')
            print(' ')

    def print_key_info(self,keys):
        for f in self.files:
            for key in keys:
                print(f[key],end='  ')
            print(' ')

    def return_full_info_as_table(self):
        T = Table()
        f = self.files[0];
        Nfiles = len(self.files)
        for key in f.keys():
            keyArray = np.empty(Nfiles,dtype = 'S100')
            for i,f2 in enumerate(self.files):
                keyArray[i] = f2[key]
            T[key] = keyArray
        return T        
                

################################################################################

def download_fitsfile_from_gemini_archive(f, downloadDir, overwrite=0, verbose=0, cookie = None):
    downloadUrl = "https://archive.gemini.edu/file/"
    if '.bz2' in f:
        myf = f.split('.bz2')[0]
    else:
        myf = f 
    durl = os.path.join(downloadUrl,myf)
    #download file
    targetFileName = os.path.join(downloadDir,myf)
    if ( (not os.path.isfile(targetFileName)) or (overwrite == 1) ):
        print('Dowloading %s'%  myf)
        testfile = urllib.URLopener()
        if cookie is not None:
            testfile.addheader('cookie', '%s=%s'% (cookie['name'],cookie['value']))
        testfile.retrieve(durl,targetFileName)
    elif verbose:
        print('File already downloaded: %s' % myf)        

def select_matching_calibration_from_gemini_archive(baseUrl,basequeryString,tdate):
    
    timeWindowSize = 0.*u.day       
    t1,t2 = get_time_window(tdate,timeWindowSize)

    Nfiles = 0
    
    # adding time constraint on selection     
    while Nfiles==0:
        
        timeString = '/%s-%s' % (t1,t2)
        queryString = basequeryString + timeString
        url = baseUrl + queryString
        gsum = JsonSummary(url)
    
        Nfiles = len(gsum.files)
        if Nfiles == 0:
            timeWindowSize += 1.*u.day  
            t1,t2 = get_time_window(tdate,timeWindowSize)
            
            #       if no files are found, increase time window around obs date
            print('No files found in time window, repeat with larger one')

    Tg = gsum.return_full_info_as_table()
    return Tg

def make_gemini_master_bias(Tg, dataDir, overwrite=0, verbose=0):
    masterbiasFitsFile = os.path.join(dataDir,"bias.fits")
    if ( (not os.path.isfile(masterbiasFitsFile)) or (overwrite == 1) ):
        Tg[['filename']].write(os.path.join(dataDir,'bias.lis'),format='ascii.no_header',formats={'filename': lambda s: s[0:-4]})                    
        gmos.gbias("@%s" % (os.path.join(dataDir,'bias.lis')), outbias=masterbiasFitsFile, rawpath=dataDir, fl_vardq=iraf.yes) 
        
        # clean up
        for fname in Tg['filename']:
            homeDir = os.environ['HOME']
            os.remove(os.path.join(homeDir,'g'+fname[0:-4]))
        os.remove(os.path.join(homeDir,'gmos.log'))
        
    elif verbose:
        print('Master BIAS already exists: %s' % masterbiasFitsFile)
        
    return masterbiasFitsFile   

def make_gemini_master_flat(Tg, dataDir, masterbiasFitsFile, overwrite=0, verbose=0):
    masterflatFitsFile = os.path.join(dataDir,"flat.fits")
    if ( (not os.path.isfile(masterflatFitsFile)) or (overwrite == 1) ):
        Tg[['filename']].write(os.path.join(dataDir,'flat.lis'),format='ascii.no_header',formats={'filename': lambda s: s[0:-4]})                    
        gmos.giflat("@%s" % (os.path.join(dataDir,'flat.lis')), outflat=masterflatFitsFile, rawpath=dataDir+'/', fl_vardq=iraf.yes, bias=masterbiasFitsFile, combine="median") 

        # clean up
        for fname in Tg['filename']:
            homeDir = os.environ['HOME']
            os.remove(os.path.join(homeDir,'g'+fname[0:-4]))
            os.remove(os.path.join(homeDir,'rg'+fname[0:-4]))
        os.remove(os.path.join(homeDir,'gmos.log'))
    elif verbose:
        print('Master FLAT already exists: %s' % masterflatFitsFile)
        
    return masterflatFitsFile

def make_gemini_frame_mosaic_PDF(filename, overwrite=0,outputDir=None, deleteFits=True):
    dataDir, fname = os.path.split(filename)
    if outputDir is not None:
        dataDir = outputDir     
    mosaicName = os.path.join(dataDir,'m'+fname)    
    saveFile = mosaicName.replace('.fits','.pdf');
    if ( (not os.path.isfile(saveFile)) or (overwrite == 1) ):          
        gmos.gmosaic(filename,outimages=mosaicName)
        fig = pl.figure(figsize=(7, 7),facecolor='w', edgecolor='k'); pl.clf();        
        gc = aplpy.FITSFigure(mosaicName, figure=fig)#,convention='wells')#, dimensions=[1 ,0]);
        gc.show_grayscale(invert = False)#, stretch='log', vmid=-1)#,vmid=-1)#, aspect = pixScaleAC_mas/pixScaleAL_mas, pmax =90 )
        if (('mbias.pdf' in saveFile) | ('mflat.pdf' in saveFile)):
            gc.hide_ytick_labels()
            gc.hide_xtick_labels()      
        gc.save(saveFile, dpi=300);
        if deleteFits:
            os.remove(mosaicName)

def make_gemini_frame_mosaic(filename, overwrite=0, outputDir = None, deleteFits = False, makePDF = False, outputDirPDF = None):
    dataDir, fname = os.path.split(filename)
    if outputDir is not None:
        dataDir = outputDir     
    mosaicName = os.path.join(dataDir,'m'+fname)    
    
    if ( (not os.path.isfile(mosaicName)) or (overwrite == 1) ):
        gmos.gmosaic(filename, outimages=mosaicName)

        # 2017-12-12 do not interpolate when creating the mosoaic
        # see http://www.gemini.edu/sciops/data/IRAFdoc/gmosinfo.html
        # gmos.gmosaic(filename, outimages=mosaicName, geointer="nearest")


    if makePDF:
        if outputDirPDF is not None:            
            saveFile = os.path.join(outputDirPDF,'m'+fname).replace('.fits','.pdf'); 
        else:
            saveFile = mosaicName.replace('.fits','.pdf');
        if ( (not os.path.isfile(saveFile)) or (overwrite == 1) ):          
            fig = pl.figure(figsize=(7, 7),facecolor='w', edgecolor='k'); pl.clf();        
            gc = aplpy.FITSFigure(mosaicName, figure=fig)#,convention='wells')#, dimensions=[1 ,0]);
            gc.show_grayscale(invert = False)#, stretch='log', vmid=-1)#,vmid=-1)#, aspect = pixScaleAC_mas/pixScaleAL_mas, pmax =90 )
            if (('mbias.pdf' in saveFile) | ('mflat.pdf' in saveFile)):
                gc.hide_ytick_labels()
                gc.hide_xtick_labels()      
            gc.save(saveFile, dpi=300);
            
    if deleteFits:
        os.remove(mosaicName)