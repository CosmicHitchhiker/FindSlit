# find_slit
GUI for determining the precise postion of long-slit spectrum.

## Purpose and limits
In work with long-slit spectrosopy data it is crucial to know the position of each spectral "row" 
(1-D spectrum that corresponds to the specific postion at the slit). It is not always possible to
determine it during observations due to uncertancies in the telescope positioning, abberations in the
otical system of a spectrograph, later correstions in the reduction etc. So, one of the important steps
of the spectral data processing is to perform the astrometrical calibration of the observartion results.
The goal of this software is to make it possible to perform this calibration in a user-friendly way.

This program uses photometry with 

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

