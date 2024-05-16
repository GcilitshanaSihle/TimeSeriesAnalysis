import glob
import numpy as np
from astropy.io import fits
from astropy.time import Time
import os 
import os.path 
import math 
import pickle 
import numpy as np 
from numpy import fft 
import glob
import pylab
import random
import re
import sys
import itertools
from astropy.coordinates import SkyCoord 
from astropy.io import fits 
from astropy.wcs import WCS
from astropy.time import Time 
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import date
# Add VAriable sources
# add on source 230 = (2615.448 px, 1908.028 px)
# source 1731 = 4677.948 px, 3970.529 px

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def give_edges(locations, nx, ny, nx_psf, ny_psf):

    nx0 = nx_psf//2
    ny0 = ny_psf//2

    edges = []

    for location in locations:
        p, q = location

        # Calculate edges for each location
        dxl = p - nx0
        xl = np.maximum(dxl, 0)

        dxu = p + nx0
        xu = np.minimum(dxu, nx)

        dyl = q - ny0
        yl = np.maximum(dyl, 0)

        dyu = q + ny0
        yu = np.minimum(dyu, ny)

        xlpsf = np.maximum(nx0 - p , 0)
        xupsf = np.minimum(nx0 + nx - p, nx_psf)

        ylpsf = np.maximum(ny0 - q, 0)
        yupsf = np.minimum(ny0 + ny - q, ny_psf)

        edges.append((slice(xl, xu), slice(yl, yu), \
                      slice(xlpsf, xupsf), slice(ylpsf, yupsf)))
        
    edges =np.array(edges)
    return edges

def gaussian(x, amp, mean, std_dev):
    return amp * np.exp(-((x - mean)**2) / (2 * std_dev**2))


def process_images():
    
    Residual_fits = glob.glob("/home/gcilitshana/Project/New/madroon-1/obs-u0-qc2/im2htc-cad1/images_original/cad1-t*-image.fits")
    PSF_fits =  glob.glob("/home/gcilitshana/Project/New/madroon-1/obs-u0-qc2/im2htc-cad1/images_original/cad1-t*-psf.fits")
    Residual_fits = natural_sort(Residual_fits)
    PSF_fits = natural_sort(PSF_fits)

    timestamps = []
    for images in (Residual_fits):
        ff = fits.open(images)
        hdr = ff[0].header
        map_date = hdr.get('DATE-OBS')
        t_mjd = Time(map_date, format='isot', scale='utc').mjd
        tmjd = round(float(t_mjd))
        timestamps.append(t_mjd)

    timestamps = np.array(timestamps)

    # Generate 10 Gaussian profiles
    np.random.seed(42)
    num_profiles = 10  # Number of Gaussian profiles
    # Initialize an empty list to store the profiles
    gaussian_profiles = []
    for _ in range(num_profiles):
        t = np.linspace(0, 3000, 3000)
        # Generate a Gaussian profile with 9 peaks
        num_peaks = 7
        amp = np.random.uniform(0.0005,0.0008, num_peaks)  # Random amplitudes
        mean = np.linspace(2156,2169, num_peaks)    # Equally spaced means between 100 and 2400
        std_dev = np.random.uniform(0.5, 1.5, num_peaks)  # Random standard deviations between 10 and 100
        profile = np.sum([gaussian(t, amp[i], mean[i], std_dev[i]) for i in range(num_peaks)], axis=0)
        gaussian_profiles.append(profile)

    gaussian_profiles =np.array(gaussian_profiles)
    
    locations = [(1853,2665),(2615, 4094), (2557, 3140),(1389,3926),(3319,3158)]  # Add more locations if needed
    nx =  6000
    ny =  6000
    nx_psf = 6000
    ny_psf = 6000
    ntime = 2942

    edges = give_edges(locations, nx, ny, nx_psf, ny_psf)

    for (i, (res, psf)) in enumerate(zip(Residual_fits, PSF_fits)):
        ff1 = fits.open(psf)
        hdr1 = ff1[0].header
        img1 = ff1[0].data[0, 0][::-1, :].T
        ff = fits.open(res)
        hdr = ff[0].header
        img = ff[0].data[0, 0][::-1, :].T
        
        for j in range(5):
            Ix, Iy, Ix_psf, Iy_psf = edges[j]
            img[Ix, Iy] += img1[Ix_psf, Iy_psf] * gaussian_profiles[j,i]

        hdu = fits.PrimaryHDU(header=hdr1)
        img = img.T[::-1, :][None, None]
        img = np.require(img, dtype=img1.dtype, requirements='F')
        hdu.data = img
        hdu.writeto(f'/home/gcilitshana/Project/New/madroon-1/obs-u0-qc2/im2htc-cad1/images/cad1-t{str(i).zfill(4)}-image.fits', overwrite=True)
process_images()