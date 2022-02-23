import app.imgscale
from astropy.io import fits
from PIL import Image
from termcolor import colored


def convertImage(file, fileName, contrast, saveDir):
    image_data = fits.getdata(file)
    jpg_file_name = fileName + '.jpg'
    if len(image_data.shape) == 2:
        sum_image = image_data
    else:
        sum_image = image_data[0] - image_data[0]
        for single_image_data in image_data:
            sum_image += single_image_data

    # sum_image = lib.imgscale.sqrt(sum_image, scale_min=0, scale_max=np.amax(image_data))
    sum_image = app.imgscale.sqrt(sum_image)
    sum_image = sum_image * contrast
    im = Image.fromarray(sum_image)
    if im.mode != 'RGB':
        im = im.convert('RGB')

    im.save(saveDir + '/' + jpg_file_name)
    im.close()
    print('[', colored('OK', 'green'), '] Converted image (', jpg_file_name.strip(), ') save')
