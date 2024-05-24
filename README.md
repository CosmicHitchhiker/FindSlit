# find_slit
GUI for determining the precise postion of a long-slit spectrum.

## Purpose of this program
In work with long-slit spectrosopy data it is crucial to know the position of each spectral "row" 
(1-D spectrum that corresponds to the specific postion at the slit). It is not always possible to
determine it during observations due to uncertancies in the telescope positioning, abberations in the
otical system of a spectrograph, later correstions in the reduction etc. So, one of the important steps
of the spectral data processing is to perform the astrometrical calibration of the observartion results.
The goal of this software is to make it possible to perform this calibration in a user-friendly way.

This program works with an image of an object or a sky region (photometry) and a spectrum of the same object.
Both files should be in the fits fomat and there are several requerments to them that are listed below.
The program finds the position of the slit where the total flux measured by spectroscopy (summ of the fluxes
at each wavelength) has maximum correlation to the flux at the image at this slit position. It is also possible
not to use all of the wavelengths but to specify the range, this can be useful when the image is observed with
some other instrument.

## Requerements to the files

#### Spectrum
The spectrum of and object should be presented as a fits-file with a 2D spectrum in its first HDU. It is considered
that the positional coordinate along the slit is the Y-coordinate (vertical axis, has index 2 in the header keywords),
and the wavelength coordinate is the X-coordinate (horizontal, has index 1 in the keywords). The flux at each pixel
can be set in any linear units (e.g. ADU/pix, erg/s/cm^2, but not the magnitude).

The necessary keywords for the header of this HDU are:
 - RA
 - DEC
 - NAXIS2

The optional, but important keywords are:
 - EPOCH
 - CDELT2
 - CRPIX2
 - POSANG
 - PARANGLE and ROTANGLE
 - SLIT
 - CDELT1, CRPIX1, CRVAL1

#### Photometry
The photometry of an object/sky region also shoud be a fits-file with a 2D image in its first HDU. It should
have correct WCS keywords in the header. One of the best ways to obtain fits-file with correct WCS keywords is to
use astrometry.net for the astrometrical calibration.

## Installation
Currently, this software was only tested on Linux systems and Python 3.10.

First of all, you need to install the necessary packages:
 - numpy
 - matplotlib
 - scipy
 - astropy
 - reproject
 - PySide6

You can do it in any way you prefer. The easiest way (but not the best one) is this command:

`pip install numpy matplotlib scipy astropy reproject PySide6`

Then you can clone this repository in any folder and make the main script executable:

```
git clone https://github.com/CosmicHitchhiker/FindSlit.git
cd FindSlit
chmod a+x find_slit.py
```
## Usage

