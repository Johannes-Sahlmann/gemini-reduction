"""
Reduce GMOS-S and GMOS-N images

This script implements:
- download data from Gemini archive given an object name
- associate and download corresponding calibrations (bias, flat)
- generate the master calibration files
- reduce the imaging data using the master calibrations
- generate pdf images

J. Sahlmann, STScI/AURA, 2016-12-12

Programmatic interface to Gemini archive implemented with the help of
https://archive.gemini.edu/help/api.html

Requires pyraf/iraf




"""

# Load the required packages
from __future__ import print_function
import os, sys
import pickle
import numpy as np
from astropy.table import Table, Column
import pylab as pl
from astropy.time import Time, TimeDelta
import glob
import shutil

from pyraf import iraf
from iraf import gemini
from pyraf.iraf import gmos

iraf.unlearn(gemini)
iraf.unlearn(gmos)




import gemini_reduction as gred
reload(gred)




# if data not public, read cookie information from the .netrc file in your home directory,
data_is_public = True
if data_is_public is False:
    import netrc
    HOST = 'gemini'
    secrets = netrc.netrc()
    username, account, cookie = secrets.authenticators( HOST )

    cookie = {'name': 'gemini_archive_session', 'value': cookie}
else:
    cookie = None


home_dir = os.environ['HOME']


demo_only = True


overwrite = 0
overwrite_archive_file_summary = 1
overwrite_calibration_association = 1
overwrite_archive_file_summary = 0
overwrite_calibration_association = 0

base_dir = os.path.join(home_dir, 'gemini_images')

object_header_name = '2MASS 1059-2113'
identifier = object_header_name.replace(' ', '')
instrument = 'GMOS-S'

makePDF = True
# makePDF = False  # aplpy version trouble 2019-06-17

################################################################################

# directory to write file summaries
tmp_dir = os.path.join(base_dir, instrument, 'results', identifier)

dataDir = os.path.join(base_dir, instrument, 'data')
rawDir = os.path.join(dataDir,'raw')
calibDir = os.path.join(dataDir,'calibrations')
reducedDir = os.path.join(dataDir,'reduced', identifier)
mosaicDir = os.path.join(reducedDir,'mosaics')
thumbnail_dir = os.path.join(reducedDir, 'thumbnails')
gred.make_dir(calibDir)
gred.make_dir(reducedDir)
gred.make_dir(mosaicDir)
gred.make_dir(thumbnail_dir)
gred.make_dir(tmp_dir)


################################################################################
# CHECK that raw data is up to date with archive, xmatch data on disk with jsonsummary

# from: https://archive.gemini.edu/help/api.html
# Construct the URL. We'll use the jsonfilelist service
baseUrl = "https://archive.gemini.edu/jsonsummary/canonical/"

################################################################################


# GET FILE SUMMARY FROM ARCHIVE

if instrument == 'GMOS-N':
    queryString = '%s/object=%s' % (instrument,object_header_name.replace('+','%252B').replace(' ','+'))
else:
    queryString = '%s/object=%s' % (instrument,object_header_name.replace(' ','+'))

queryString += '/OBJECT'
queryString += '/science'
queryString += '/imaging'
queryString += '/Win' # added 2017-11-01


retrieve_only_recent_data = False
if retrieve_only_recent_data:
    # add constraint on time to avoid 250 entry limit of the Gemini archive interface
    from datetime import datetime
    t_now = Time(str(datetime.now()), format='iso')
    t_m1yr = t_now - TimeDelta(365, format='jd')
    queryString += '/%s-%s'%(t_m1yr.iso.split()[0].replace('-',''), t_now.iso.split()[0].replace('-',''))


url = str(baseUrl + queryString)

geminiArchiveFileSummaryFile = os.path.join(tmp_dir, '%s_geminiArchiveFileSummary.txt' % identifier)
if ( (not os.path.isfile(geminiArchiveFileSummaryFile)) or (overwrite_archive_file_summary == 1) ):
    ga = gred.JsonSummary(url)
    T = ga.return_full_info_as_table()
    T.write(geminiArchiveFileSummaryFile,format='ascii.basic', overwrite=True)
else:
    T = Table.read(geminiArchiveFileSummaryFile,format='ascii.basic')

# identify number of observing blocks
observation_ids = np.unique(T['observation_id'])
NOB = len(observation_ids)
print('*'*60)
print('*'*60)

print('%s has %d exposures in Gemini Archive in %d observing blocks' % (identifier, len(T), NOB))


for f in T['filename']:
    gred.download_fitsfile_from_gemini_archive(f, rawDir, overwrite, verbose=1, cookie=cookie)

    if demo_only:
        # stop after first file for demo purposes
        break
    
geminiArchiveFileSummaryFile_calibAssoc = os.path.join(tmp_dir, '%s_geminiArchiveFileSummary_calibAssoc.txt' % identifier)

if ( (not os.path.isfile(geminiArchiveFileSummaryFile_calibAssoc)) or (overwrite_calibration_association == 1) ):

    T['masterbiasFitsFile'] = np.array(['None']*len(T)).astype('S300')
    T['masterflatFitsFile'] = np.array(['None']*len(T)).astype('S300')

    # main loop
    for ob_id in observation_ids:
        index_ob_id = np.where(T['observation_id']==ob_id)[0]
        baseFileIndex = index_ob_id[0]
        baseRawScienceFrameName = T['filename'][baseFileIndex][0:-4]
    
        date = T['ut_datetime'][baseFileIndex].split(' ')[0].replace('-','')
        tdate = Time(T['ut_datetime'][baseFileIndex].split(' ')[0]) 
        print('%s\t%s'%(object_header_name,ob_id))

        try:
            sd1 = T['detector_config'][baseFileIndex].split(' ')[0]
        except KeyError:
            sd1 = T['detector_gain_setting'][baseFileIndex] # e.g. "low"
        sd2 = T['detector_binning'][baseFileIndex]
        sd3 = T['detector_roi_setting'][baseFileIndex].replace(' ','')
        detectorParameterString = '/%s/%s/%s' % (sd1,sd2,sd3)
    
        calTypes = np.array(['BIAS','FLAT'])
        if instrument == 'GMOS-N':
            qa_stateString = '/NotFail'        
        else:
            qa_stateString = '/Pass'
        
        reductionStateString = '/RAW'
        observationTypeString = np.array(['/BIAS','/Twilight'])
        filterString = np.array(['','/filter=%s' % (T['filter_name'][baseFileIndex])])      
        

        masterbiasFitsFile = []                

        for j,calType in enumerate(calTypes):

            print('Looking for %s files' % calType)
            queryString = '%s' % (instrument)
            queryString += observationTypeString[j] + qa_stateString + detectorParameterString + filterString[j] + reductionStateString

            Tg = gred.select_matching_calibration_from_gemini_archive(baseUrl, queryString, tdate)

            calibration_date_root = Tg['filename'][0][1:9]

            downloadDir = os.path.join(calibDir, calibration_date_root, calType)
            gred.make_dir(downloadDir)

            for f in Tg['filename']:
                gred.download_fitsfile_from_gemini_archive(f, downloadDir, overwrite, verbose=1)
                
            if calType =='BIAS':   
                masterbiasFitsFile = gred.make_gemini_master_bias(Tg, downloadDir, overwrite, verbose=1)
    
            elif calType =='FLAT':         
                masterflatFitsFile = gred.make_gemini_master_flat(Tg, downloadDir, masterbiasFitsFile, overwrite, verbose=1)

        if makePDF:
            gred.make_gemini_frame_mosaic_PDF(masterbiasFitsFile, overwrite)
            gred.make_gemini_frame_mosaic_PDF(masterflatFitsFile, overwrite)

        T['masterbiasFitsFile'][index_ob_id] = np.array([masterbiasFitsFile]*len(index_ob_id))
        T['masterflatFitsFile'][index_ob_id] = np.array([masterflatFitsFile]*len(index_ob_id))

        if demo_only:
            break
    T.write(geminiArchiveFileSummaryFile_calibAssoc, format='ascii.basic', overwrite=True)
else:
    T = Table.read(geminiArchiveFileSummaryFile_calibAssoc, format='ascii.basic')



# reduce the science images
for j,fname in enumerate(T['filename']):
    pl.close('all')
    f = fname[0:-4]
    inFileName  = os.path.join(rawDir,f)
    outFileName = os.path.join(reducedDir,f.replace('.fits','_red.fits'))
    mosaicFileName = os.path.join(mosaicDir,'m'+f.replace('.fits','_red.fits'))
    if ( (not os.path.isfile(mosaicFileName)) or (1) ):
        if 1:
            gmos.gireduce(inimages=inFileName, outimages=outFileName, bias=T['masterbiasFitsFile'][j], flat1=T['masterflatFitsFile'][j])

            # clean up
            homeDir = os.environ['HOME']
            os.remove(os.path.join(homeDir,'g'+f))
            os.remove(os.path.join(homeDir,'gmos.log'))

        gred.make_gemini_frame_mosaic(outFileName, overwrite=0, outputDir = mosaicDir, makePDF=makePDF, outputDirPDF = thumbnail_dir)
        # clean up
        os.remove(outFileName)
    else:
        print('Mosaic already exists: %s' % mosaicFileName)

    if demo_only:
        break

sys.exit()


