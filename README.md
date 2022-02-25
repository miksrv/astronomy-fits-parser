# Astronomical Image Analyzer 

This script is designed to analyze astronomical data files (**FITS**). 
The algorithms used were written by several authors, this script combines all the tools into one and provides the following features:
- Calculation of the background of the starry sky and search for stars in the image.
- Computing the FHWM for each star in the image and finding the average.
- Calculation of the signal-to-noise parameter for each star and the average value of the entire field.
- Creation of report files (txt and csv) for each FITS file.
- Plotting Calculated Values.
- Convert FITS to JPEG image file and send to remote server.
- Sending FITS file headers and forced values (number of stars, FWHM, SNR, etc.) to a remote server.

![GitHub Light](docs/script-running.png)

This script is used to analyze captured frames of astronomical images from a homemade astronomical observatory. 
File analysis data are uploaded to a remote server, where they are visualized in order to accumulate and store information about the astronomical objects being photographed. 
Using the calculated parameters (FWHM, etc.), you can evaluate the quality of captured frames, and image preview allows you to evaluate frames immediately after uploading to the server.

An example of uploaded data can be viewed here: 
https://observatory.miksoft.pro/object/IC_1396

The repository for this service is here: https://github.com/miksrv/observatory

### Script authors
1. [Teresa Symons](https://github.com/tasymons/photparty)
2. [Denis Bernardes](https://github.com/DBernardes/)
3. [김동현](https://github.com/psds075)

### Need help!
Please help improve this script. There are two important tasks, well, maybe at your professional discretion:
1. Optimize the application cycle.
2. Optimize the calculation of the FWHM of stars - for example, transfer a data array with the area around the star, 
because now the source file is opened every time for FWHM analysis, which is very resource intensive. 

## Installing and running the script 
The script is written in Python (v 3.10), to run it, you need to have a language interpreter installed on your computer. 
On an OS such as Linux, it comes by default, and for Windows, you need to download and install from the official site: 
https://www.python.org/downloads/windows/

1. Go to the directory with the downloaded script, open console and install necessary additional packages. 
Type in the command line:
`pip install numpy configparser matplotlib astropy scipy`
2. Set up a configuration file. Set the following values in the **config.ini** file:
   - Path to the root directory with FITS files. The script works recursively and checks all nested directories:
`path = /home/user/fits`
   - If the FWHM calculation option is enabled `calculateFWHM = on`, then you need to specify the path to the bias calibration file: 
`calibrateBiasFile = /home/user/calibrate/bias.fits`
   - If a remote service is configured to receive data about analyzed files `toAPI = on`, then you need to specify the API:
`toAPIEndpoint = http://api.miksoft.pro/astro/set/fits`
   - If you want to upload image files that will be converted from FITS `upload = on`, then you must also specify the API:
`uploadAPI = http://api.miksoft.pro/astro/set/image`
3. The minimum configuration is done, you need to run the script through the console (command line):
` py.exe .\main.py`
4. After parsing each file, if the corresponding options are enabled in the settings:
   - The image converted from the FITS file JPG will be added to the directory (`/images/`).
   - The `/analysis/` directory will contain two files for each FITS file.
   - A `.lastdate` file will be created at the root of the program, which contains the date of the last processed file. Thus, running the script again will not analyze already processed FITS files.
   - Image data, titles, stars, etc. will be sent to the API.

## Configuration Options
In the root directory with the program there is a configuration file `config.ini`, below are a description of its parameters.

### [GENERAL]
Specify the path to the folder containing the files for analysis (FITS). All contained directories will also be scanned. 

`path = /home/user/fits`

Calculate the average FWHM value for all found stars in the frame 
(this may take time, depending on the power of the computer and the number of stars found in the frame). 
Type 'on' or 'off'

`calculateFWHM = on`

Calibration frame (bias) containing dark current noise. Used to calculate the SNR of a star.
A bias image is an image acquired with the same CCD operation mode used to acquire the star image, with no
light incidence, and for an exposure time equal to zero.

`calibrateBiasFile = /home/user/calibrate/bias.fits`

### [REPORT]
If this option is enabled, a report will be generated for each parsed file in the form of two text files.
The first file contains the coordinates of the detected stars and the calculated brightness,
the second file contains detailed information on the analysis. Type 'on' or 'off'

`toFiles = on`

Directory for saving text reports. Can be either absolute or relative (from the directory with the script)

`toFilesPath = analysis`

Send a report to a remote server via API or not.
FITS file header + parameters for calculating the number of stars + FHWM averages will be sent (if calculation is enabled)

`toAPI = on`

Server API address. All data is sent by POST request as a JSON array.
In case of success the server should return a response with a 200 header and JSON data with a (boolean) parameter { status: true }

`toAPIEndpoint = http://api.miksoft.pro/astro/set/fit_object`

### [IMAGE]
Convert FITS file to image JPG? Type 'on' or 'off'

`convertJPG = on`

The contrast (brightness) value of the image. The higher the number,
the brighter the image will be with more detail. It is recommended to set 300-500

`contrast = 500`

Directory for saving converted images. Can be either absolute or relative (from the directory with the script)

`saveDir = images`

If the option to convert images (FITS to JPG) is enabled, then this option will allow uploading to a remote server

`upload = on`

Remote API for upload JPG image, should be return (if success) 200 header and JSON data with a parameter { status: true }

`uploadAPI = http://api.miksoft.pro/astro/set/photo`

### [PARSING]

Parameters for background sampling: width/length in pixels of box for random background sampling to determine background value

`backSize = 5`

Number of random background samples to take

`backNum = 1000`

Selection of area of each frame to analyze:
 - If area of interest is in the central 50% of frame, select 'half'
 - If area of interest is entire frame, select 'whole'
 - If custom area of interest is desired, select 'custom' and define range of X and Y coordinates with (1,1) as the bottom left corner of the frame
 
```
frameArea = whole
xLow = 900
xHigh = 1500
yLow = 1200
yHigh = 1750
```

Pixel rejection: select pixel value above which values will be zeroed out

`upLim = 5000`

Select number of sigma below which negative pixel values will be zeroed out. 
Example: lowsig = 3 means pixel values less than -3*inset standard deviation will become 0

`lowSig = 3`

Detection level: select number of standard deviations above background required for star detection

`sig = 20`

Square-aperture size: select the half-width of the box used for photometry

`boxHW = 15`

### [HEADER]

Header key words/values. if no keyword exists, enter 'NONE' and define value instead.

Exposure time, Filter, Airmass, Gain
```
exptimeKey = EXPTIME
exptimeVal = 60
filterKey = FILTER
filterVal = b
airmassKey = AIRMASS
airmassVal = 1.0
gainKey = GAIN
gainVal = 90
```

### [PLOTS]

Suppress or display plots of summed rows/columns with detection level marked (on or off)

`plotDetect = 'off'`

Suppress or display plots of fits image with detected star apertures overlaid (on or off)

`plotStars = 'off'`