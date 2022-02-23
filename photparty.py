# The goal of this script is to read fits files, locate stars, and return magnitudes for those stars
# Background sky level is calculated using random median sampling of square arrays
# A specific inset of the data array is taken if desired
# Pixel values are summed row-wise and column-wise and background level is used to determine which rows and columns contain stars
# A list of matched star coordinates is prepared
# Star magnitudes are obtained using square apertures after being background subtracted for the region specific to the star
# By Teresa Symons 2016


import numpy as np
import os
import configparser
import matplotlib.patches as patches
import matplotlib.pylab as plt
from pathlib import Path
from termcolor import colored
from astropy.io import fits
from astropy.table import Table
from app import FWHM
from app.progressBar import progressBar
from app.convertImage import convertImage
from app.binsum import binsum
from app.background import background
from app.starlocate import starlocate
from app.starmed import starmed
from app.starphot import starphot

config = configparser.ConfigParser()
config.read('config.ini')

# Читаем дату последнего обработанного файла
lastDate = 0
lastDateFile = '.lastdate'

# Если файл с последней датой существует, читаем его и переписываем переменную с датой
if os.path.exists(lastDateFile):
    tmpFile = open(lastDateFile)
    tmpStr = tmpFile.read().strip()
    lastDate = 0 if not tmpStr else float(tmpStr)

# Рекурсивно ищем файлы с определенным расширением
filesList = [path for path in Path(config['GENERAL']['path']).rglob('*.fits') if lastDate < os.stat(path).st_ctime]
filesList = sorted(filesList, key=os.path.getctime)
totalFiles = len(filesList)

if totalFiles == 0:
    print(colored('No new files found in work directory', 'green'))
    quit()

# Количество обработанных файлов
iterationCount = 0

for file in filesList:
    iterationCount += 1

    # --- СТАРТ ЦИКЛА --- #
    fileName = file.name
    filePath = os.path.dirname(file)

    print(colored('[ ' + str(iterationCount) + '/' + str(totalFiles) + ' ] Start parsing FITS: ' + fileName, 'blue', None, ['bold']))

    # Create and open analysis files if reports enabled in config
    if config['GENERAL']['report'] == 'on':
        fNameStar = config['GENERAL']['reportPath'] + '/' + fileName + '.stars.csv'
        fNameData = config['GENERAL']['reportPath'] + '/' + fileName + '.data.txt'
        fileStar = open(fNameStar, 'w')
        fileData = open(fNameData, 'w')

    # FITS Open file
    image = fits.open(file)

    # Retrieve data, exposure time, airmass, and filter from fits file:
    Data = image[0].data

    etime = config['HEADER']['exptimeVal'] if config['HEADER']['exptimeKey'] == 'NONE' else image[0].header[config['HEADER']['exptimeKey']]
    filter = config['HEADER']['filterVal'] if config['HEADER']['filterKey'] == 'NONE' else image[0].header[config['HEADER']['filterKey']]
    airmass = config['HEADER']['airmassVal'] if config['HEADER']['airmassKey'] == 'NONE' else image[0].header[config['HEADER']['airmassKey']]
    gain = config['HEADER']['gainVal'] if config['HEADER']['gainKey'] == 'NONE' else image[0].header[config['HEADER']['gainKey']]

    # Compute background sky level through random median sampling:
    # Inputs: data array, nxn size of random subarray to use for sampling, and number of desired sampling iterations
    back, skyvals = background(Data, int(config['PARSING']['backSize']), int(config['PARSING']['backNum']))

    # Create desired inset of total data array:
    if config['PARSING']['frameArea'] == 'half':
        # Find midpoint of image data
        mid = len(Data) / 2
        # Take an inset of data that is half of area
        inset = Data[round(mid / 2):3 * round(mid / 2), round(mid / 2):3 * round(mid / 2)]
        xlow = 0
        xhigh = 0
        ylow = 0
        yhigh = 0
    if config['PARSING']['frameArea'] == 'whole':
        # Use entire data array as inset
        inset = Data
        mid = 0
        xlow = 0
        xhigh = 0
        ylow = 0
        yhigh = 0
    if config['PARSING']['frameArea'] == 'custom':
        inset = Data[
            int(config['PARSING']['yLow']) - 1:int(config['PARSING']['yHigh']) - 1,
            int(config['PARSING']['xLow']) - 1:int(config['PARSING']['xHigh']) - 1
        ]
        mid = 0

    # Blanket removal of bad pixels above 45000 and 3*standard deviation below 0:
    inset[inset > int(config['PARSING']['upLim'])] = 0
    std = np.std(inset)
    inset[inset < -int(config['PARSING']['lowSig']) * std] = 0

    # Calculate sky background for specific inset:
    # Inputs: inset data array, nxn size of random subarray used for sampling, number of desired sampling iterations
    insetback, insetskyvals = background(inset, int(config['PARSING']['backSize']), int(config['PARSING']['backNum']))

    # Compute summed row and column values for desired array by number of bins:
    # Inputs: inset data array, number of bins desired
    rowsum, colsum = binsum(inset, 1)

    # Locate values in summed row and column vectors that are greater than desired sigma level above background:
    # Inputs: Data array, background level variable, desired sigma detection level, summed row vector, summed column vector
    starrow, starcol, backsum, std, sigma = starlocate(inset, insetback, int(config['PARSING']['sig']), rowsum, colsum)
    if not starrow:
        print('[ ' + colored('ERROR', 'red'), '] No stars found in ' + fileName + ' by row - check threshold.')

    if not starcol:
        print('[ ' + colored('ERROR', 'red'), '] No stars found in ' + fileName + ' by column - check threshold.')

    if starrow != [] and starcol != []:

        # Plot summed row and column values with detection level marked:
        if config['PLOTS']['plotDetect'] == 'on':
            plt.plot(rowsum)
            plt.plot((0, len(rowsum)), (backsum + sigma * std, backsum + sigma * std))
            plt.title('Summed Rows' + '-' + fileName)
            plt.xlabel('Row Index in Data Inset')
            plt.ylabel('Summed Row Value')
            plt.show()
            plt.plot(colsum)
            plt.plot((0, len(colsum)), (backsum + sigma * std, backsum + sigma * std))
            plt.title('Summed Columns' + '-' + fileName)
            plt.xlabel('Column Index in Data Inset')
            plt.ylabel('Summed Column Value')
            plt.show()

        # Take indices of detected star pixels and divide into sublists by individual star:
        # Return sublists of star indices, number of stars, and median pixel of each star
        # Pair star center row with star center column and return total number of pairs found
        # Inputs: vectors containing indices of detected star pixels for row and column and inset data array
        rowloc, colloc, numstarr, numstarc, rowmed, colmed, starpoints, adjstarpoints = starmed(
            starrow, starcol, inset, mid, xlow, ylow
        )

        # Take list of star coordinates and find summed pixel values within a square aperture of desired size:
        # Also find background values for each star and subtract them from star values
        # Convert background-subtracted star values into fluxes and then magnitudes
        # Inputs: half-width of nxn square aperture, inset data array, vector containing star coordinates, exposure time
        boxsum, starback, backsub, flux, mags, hw, magerr = starphot(
            int(config['PARSING']['boxHW']), inset, starpoints, etime, gain, fileName
        )

        # Output data table to file if reports enabled in config:
        if config['GENERAL']['report'] == 'on':
            n = len(mags)
            tname = [fileName] * n
            tfilter = [filter] * n
            tairmass = [airmass] * n
            tetime = [etime] * n
            x = [x for [x, y] in adjstarpoints]
            y = [y for [x, y] in adjstarpoints]
            t = Table([tname, tfilter, tairmass, tetime, x, y, mags, magerr], names=('File_Name', 'Filter', 'Airmass', 'Exposure_Time', 'X', 'Y', 'Magnitude', 'Mag_Err'))
            t.write(fileStar, format='csv')

        # Plot fits image with square apertures for detected stars overlaid
        if config['PLOTS']['plotStars'] == 'on':
            fig, ax = plt.subplots(1)
            ax.imshow(Data, cmap='Greys', vmin=0, vmax=10)
            for i in range(0, len(x)):
                rect = patches.Rectangle(((x[i] - 1 - hw), (y[i] - 1 - hw)), 2 * hw, 2 * hw, linewidth=1, edgecolor='r', facecolor='none')
                ax.add_patch(rect)
            if config['PARSING']['frameArea'] == 'custom':
                rect = patches.Rectangle((xlow - 1, ylow - 1), xhigh - xlow, yhigh - ylow, linewidth=1, edgecolor='b', facecolor='none')
                ax.add_patch(rect)
            plt.title(fileName)
            plt.show()

        file_dir = os.path.dirname(os.getcwd()) + '\\\photparty\\bias\\'
        bias_name = 'Bias_0.000032_secs_2021-08-25T22-24-49_015.fits'

        counter = 0
        progress = 0
        meanFWHM = 0
        meanSNR = 0

        # Initial call to print 0% progress
        totalLength = len(adjstarpoints)

        if config['GENERAL']['calculateFWHM'] == 'on':
            if totalLength != 0:
                progressPrefix = '[ ' + colored('FWHM', 'yellow') + ' ] Calculation'
                progressBar(0, totalLength, prefix=progressPrefix, suffix='', length=50)

                for [x, y] in adjstarpoints:
                    try:
                        FWHM_obj = FWHM.fwhm(img_name=file, xy_star=(x, y), sky_radius=30, bias_name=file_dir + '\\' + bias_name, ccd_gain=gain)
                        FWHM_obj.read_star_img()
                        FWHM_obj.get_max_count()
                        FWHM_obj.set_centroid()
                        fwhm, star_radius, x, y = FWHM_obj.calc_FWHM()

                        FWHM_obj.read_bias_img()
                        FWHM_obj.calc_dark_current()
                        FWHM_obj.read_exp_time()
                        FWHM_obj.read_em_gain()
                        FWHM_obj.calc_star_sky_flux()
                        snr, rn, sky_flux, star_flux, n_pixels, bias_level = FWHM_obj.calc_SNR()

                        counter += 1
                        meanFWHM += fwhm
                        meanSNR += snr

                        # print('FWHM: ', round(fwhm, 2), 'pixels')
                        # print('Star Radius:', round(star_radius), 'pixels')
                        # print('Centroide: %i,%i' % (x, y))
                        # print('SNR: ', round(snr, 2))
                        # print('Sky Flux: ', round(sky_flux, 2), 'photons/s')
                        # print('Star Flux:', round(star_flux, 2), 'photons/s')
                        # print('Star Pixels:', n_pixels)

                    except:
                        continue

                    finally:
                        progress += 1
                        progressBar(progress, totalLength, prefix=progressPrefix, suffix='', length=50)

                print(
                    '[', colored('OK', 'green'),
                    '] Mean FWHM:', colored(round(meanFWHM / counter, 2), 'green', None, ['bold']),
                    '| Mean SNR:', colored(round(meanSNR / counter, 2), 'green', None, ['bold']),
                )

            else:
                print('[', colored('ERROR', 'red'), '] No stars found in file')

        print(
            '[', colored('OK', 'green'),
            '] Stars:', colored(len(starpoints), 'green', None, ['bold']),
            '| Sky background:', colored(round(insetback, 2), 'green', None, ['bold']),
            '| Deviation', colored(round(std, 2), 'green', None, ['bold'])
        )

        # Print full parsing info to file if reports is enabled in config
        if config['GENERAL']['report'] == 'on':
            print('Exposure time:', file=fileData)
            print(etime, file=fileData)
            print('Filter:', file=fileData)
            print(filter, file=fileData)
            print('Airmass:', file=fileData)
            print(airmass, file=fileData)
            print('Sky background level:', file=fileData)
            print(back, file=fileData)
            print('Shape of inset array:', file=fileData)
            print(np.shape(inset), file=fileData)
            print('Inset background sky level:', file=fileData)
            print(insetback, file=fileData)
            print('Summed background value for one row/column:', file=fileData)
            print(backsum, file=fileData)
            print('Standard deviation of inset:', file=fileData)
            print(std, file=fileData)
            print('Detection level in sigma:', file=fileData)
            print(sigma, file=fileData)
            print('Indices of detected stars:', file=fileData)
            print(starrow, file=fileData)
            print(starcol, file=fileData)
            print('Indices of detected stars divided into sublist by star:', file=fileData)
            print(rowloc, file=fileData)
            print(colloc, file=fileData)
            print('Number of stars found by row/column:', file=fileData)
            print(numstarr, file=fileData)
            print(numstarc, file=fileData)
            print('Median pixel of each star by row/column:', file=fileData)
            print(rowmed, file=fileData)
            print(colmed, file=fileData)
            print('Paired indices of star centers:', file=fileData)
            print(starpoints, file=fileData)
            print('Total number of stars found:', file=fileData)
            print(len(starpoints), file=fileData)
            print('Coordinates of stars (x,y):', file=fileData)
            print(adjstarpoints, file=fileData)
            print('Width/Height of boxes:', file=fileData)
            print(hw * 2, file=fileData)
            print('Pixel sums for boxes around each star:', file=fileData)
            print(boxsum, file=fileData)
            print('Background values for each star:', file=fileData)
            print(starback, file=fileData)
            print('Background subtracted star values:', file=fileData)
            print(backsub, file=fileData)
            print('Flux of stars:', file=fileData)
            print(flux, file=fileData)
            print('Magnitudes for each star:', file=fileData)
            print(mags, file=fileData)
            print('[', colored('OK', 'green'), '] Data files save')

    if config['IMAGE']['convertJPG'] == 'on':
        convertImage(
            file=file,
            fileName=fileName,
            contrast=int(config['IMAGE']['contrast']),
            saveDir=config['IMAGE']['saveDir']
        )
    # --- СТОП ЦИКЛА --- #

    # В САМОМ КОНЦЕ РАБОТЫ ЦИКЛА
    # Обновляем дату последнего обработанного файла
    fLastTime = open(lastDateFile, 'w')
    print(os.stat(file).st_ctime, file=fLastTime)

    print('\r')
