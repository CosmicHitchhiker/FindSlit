#! /usr/bin/env python3

# Copyright (C) 2022 The Qt Company Ltd.
# SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

import sys

import matplotlib.pyplot as plt
import numpy as np
# from scipy.stats import norm
import matplotlib
from astropy.visualization import simple_norm
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5, Angle, ICRS
import reproject
import time
from scipy import optimize
from copy import deepcopy

# from matplotlib import pyplot as plt
from matplotlib.figure import Figure
# noinspection PyUnresolvedReferences
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PySide6.QtCore import Slot, Signal, QRunnable, QThreadPool, QObject
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QDoubleSpinBox,
    QSpinBox,
    QHBoxLayout,
    QFormLayout,
    QGridLayout,
    QPushButton,
    QStatusBar,
    QFileDialog,
    QCheckBox,
    QLabel,
    QMessageBox,
    QProgressDialog,
)
from PySide6.QtGui import QKeySequence, QShortcut

from OpenFile import OpenFile
from radecSpinBox import radecSpinBox
from slitPolygon import slitPolygon, rotatedTriangle
from SetBorders import SetBorders
from SetBias import SetBias

matplotlib.use('QtAgg')


def norm_vector(vec):
    """
    Parameters
    ----------
    vec : array_like
        Input vector.

    Returns
    -------
    array_like
        Normalized vector.
    """
    try:
        return vec / np.linalg.norm(vec)
    except TypeError:
        raise TypeError('Input should be a numerical array-like object')


def correlation(vec1, vec2):
    """
    Parameters
    ----------
    vec1 : array_like
        The first vector for correlation calculation.
    vec2 : array_like
        The second vector for correlation calculation.

    Returns
    -------
    float
        The correlation between the two input vectors `vec1` and `vec2`, calculated using the arccosine
        of the dot product of the normalized vectors.

    """
    if not np.shape(vec1) == np.shape(vec2):
        raise IndexError('vectors shoud have same dimensions to find their correlation')
    a = norm_vector(vec1)
    b = norm_vector(vec2)
    return np.arccos(np.sum(a * b))


class SlitParams:
    def __init__(self):
        self.refcoord = SkyCoord(0, 0, unit=(u.hourangle, u.deg))
        self.crpix = 0.
        self.crpix_init = 0.
        self.npix = 0
        self.npix_init = 0
        self.pa = 0.
        self.scale = 1.
        # slit width (arcsec)
        self.width = 1.
        self.frame = self.refcoord.skyoffset_frame(rotation=self.pa * u.deg)
        self.pixpos = np.array([- self.crpix, self.npix - self.crpix])
        self.pos = self.pixpos * self.scale

    def from_header(self, hdr):
        """
        Parameters
        ----------
        hdr : dict
            The header containing the necessary information for calibration.

        Returns
        -------
        None

        """
        if 'EPOCH' in hdr:
            equinox = "J" + str(hdr['EPOCH'])
            # print(equinox)
            self.refcoord = SkyCoord(hdr['RA'], hdr['DEC'],
                                     unit=(u.hourangle, u.deg),
                                     frame=FK5(equinox=equinox)).transform_to(ICRS)
        else:
            self.refcoord = SkyCoord(hdr['RA'], hdr['DEC'],
                                     unit=(u.hourangle, u.deg))

        # print(self.refcoord)
        if 'CDELT2' in hdr:
            self.scale = hdr['CDELT2']  # arcsec/pix
        else:
            self.scale = 0.3

        self.npix = hdr['NAXIS2']

        if 'CRPIX2' in hdr:
            self.crpix = hdr['CRPIX2']
        else:
            self.crpix = hdr['NAXIS2'] / 2

        if 'POSANG' in hdr:
            # TDS
            self.pa = hdr['POSANG']
        elif 'PARANGLE' in hdr and 'ROTANGLE' in hdr:
            # SCORPIO
            self.pa = hdr['PARANGLE'] - hdr['ROTANGLE'] + 132.5
        else:
            self.pa = 0

        if self.scale < 0:
            self.pa += 180
            self.scale *= -1

        # if PA not in [0,360) (including negative PA)
        self.pa = self.pa - 360 * (self.pa // 360)

        if 'SLIT' in hdr:
            self.width = self.parse_tds_slit(hdr['SLIT'])
        else:
            self.width = 1.0

        self.frame = self.refcoord.skyoffset_frame(rotation=self.pa * u.deg)
        self.crpix_init = self.crpix
        self.npix_init = self.npix
        self.pixpos = np.array([- self.crpix, self.npix - self.crpix])
        self.pos = self.pixpos * self.scale

    def from_fields(self, ra, dec, pa, scale):
        self.pa = pa
        self.pos = self.pixpos * scale
        self.scale = scale
        self.refcoord = SkyCoord(ra, dec,
                                 unit=(u.hourangle, u.deg))
        self.frame = self.refcoord.skyoffset_frame(rotation=self.pa * u.deg)

    def add_limits(self, lims):
        miny, maxy = lims[2], lims[3]
        if miny is not None:
            self.crpix = self.crpix_init - miny
        else:
            self.crpix = self.crpix_init

        if maxy is not None:
            self.npix = maxy - miny + 1
        elif miny is not None:
            self.npix = self.npix_init - miny + 1
        else:
            self.npix = self.npix_init

        self.pos = np.array([- self.crpix, self.npix - self.crpix]) * self.scale

    def set_pixpos(self, pixpos):
        self.pixpos = pixpos
        self.pos = self.pixpos * self.scale

    def set_pos(self, pos):
        self.pos = pos
        self.pixpos = self.pos / self.scale

    @staticmethod
    def parse_tds_slit(slitname):
        if isinstance(slitname, (int, float)):
            return slitname

        if slitname[-4:] == 'asec':
            return float(slitname[:-4])
        # if slitname == '1asec':
        #     return 1.0
        # elif slitname == '1.5asec':
        #     return 1.5
        # elif slitname == '10asec':
        #     return 10
        print('Unknown SLIT keyword')
        return 1.0


class InterpParams:
    """Здесь повёрнутое изображение, можно брать его куски"""

    def __init__(self, image_hdu=None, slit=None):
        self.image_hdu = image_hdu
        self.bias = 0.
        if image_hdu is not None:
            self.slit_wcs = WCS(image_hdu.header)
            self.slit_hdr = image_hdu.header
            self.image_rotated = image_hdu.data
        if slit is not None:
            self.slit = slit
            self.make_slit_wcs(slit)
        self.rotate_image()

    def update(self, image_hdu=None, slit=None):
        self.image_hdu = image_hdu
        if image_hdu is not None:
            self.slit_wcs = WCS(image_hdu.header)
            self.slit_hdr = image_hdu.header
            self.image_rotated = image_hdu.data
        if slit is not None:
            self.slit = slit
            self.make_slit_wcs(slit)
        self.rotate_image()

    @staticmethod
    def get_rot_matrix(angle):
        """Matrix to transform X-Y coordinates (right-handed) to pseudo lon-lat
        (left-handed), assuming counter-clockwise rotation from Y-axis to
        lat-axis."""
        alph = np.radians(angle)
        rot_matrix = np.array([[-np.cos(alph), np.sin(alph)], [np.sin(alph), np.cos(alph)]])
        return rot_matrix

    def make_slit_wcs(self, slit, shape=None, center=None, cdelt=None):
        """Make wcs where y-axis corresponds to the slit PA.
        Refpix of the wcs is the slit center.

        Parameters
        ----------
        slit
        shape
        center
        cdelt : float
            cdelt in the resulting wcs, same for both axes, degrees

        Returns
        -------

        """
        slit_pa = slit.pa
        slitpos = slit.refcoord.fk5
        # slit length and width in arcsecs
        slit_length = slit.npix * slit.scale
        slit_width = slit.width
        slit_center_shift = slit.npix * 0.5 - slit.crpix

        if cdelt is None:
            cdelt = (0.5 * u.arcsec).to(u.deg).value
        cdeltarcsec = cdelt * 3600
        if shape is None:
            darcsec = 92
            nx = int((slit_width + 2 * darcsec)/cdeltarcsec)
            ny = int((slit_length + 2 * darcsec)/cdeltarcsec)
            shape = (nx, ny)
        if center is None:
            refpix_shift = slit_center_shift * slit.scale / cdeltarcsec
            center = [shape[0] / 2.0, shape[1] / 2.0 - refpix_shift]

        self.slit_wcs = WCS(naxis=2)
        self.slit_wcs.wcs.cdelt = [cdelt, cdelt]
        self.slit_wcs.wcs.cunit = [u.deg, u.deg]
        self.slit_wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        self.slit_wcs.wcs.pc = self.get_rot_matrix(slit_pa)
        self.slit_wcs.wcs.crpix = center
        # print('Slit center: ', slitpos)
        self.slit_wcs.wcs.crval = [slitpos.ra.to(u.deg).value, slitpos.dec.to(u.deg).value]
        # print('Slit center: ', self.slit_wcs.wcs.crval)
        self.slit_wcs.pixel_shape = shape
        self.slit_wcs.wcs.latpole = 0.
        self.slit_wcs.wcs.lonpole = 180.
        self.slit_wcs.wcs.equinox = 2000.0

        w_header = self.slit_wcs.to_header()
        w_header['NAXIS'] = 2
        w_header['NAXIS1'], w_header['NAXIS2'] = self.slit_wcs.pixel_shape

        self.slit = deepcopy(slit)
        self.slit_hdr = w_header
        return self.slit_wcs, w_header

    def rotate_image(self, slit=None):
        t = time.perf_counter()
        if slit is not None:
            # self.slit = slit
            self.make_slit_wcs(slit)
        if self.image_hdu is not None:
            self.image_rotated, _ = reproject.reproject_interp(self.image_hdu,
                                                               self.slit_hdr)
            self.image_rotated -= self.bias
            np.nan_to_num(self.image_rotated, False)
        print('reproject time ', time.perf_counter() - t)

    @Slot()
    def plot_image(self, slitpos=None):
        plt.figure()
        plt.subplot(projection=self.slit_wcs)
        norm = simple_norm(self.image_rotated, 'linear', percent=99.3)
        plt.imshow(self.image_rotated, norm=norm, cmap='bone')

        if slitpos is not None:
            high = slitpos.pos.max() * u.arcsec
            low = slitpos.pos.min() * u.arcsec
            top = SkyCoord(0 * u.deg, high,
                           frame=slitpos.frame)
            bottom = SkyCoord(0 * u.deg, low,
                              frame=slitpos.frame)
            x, y = self.slit_wcs.world_to_pixel(slitpos.refcoord)
            plt.plot(x, y, 'ro')
            bpix = [*self.slit_wcs.world_to_pixel(bottom)]
            tpix = [*self.slit_wcs.world_to_pixel(top)]
            linex, liney = [bpix[0], tpix[0]], [bpix[1], tpix[1]]

            plt.plot(linex, liney, color='tab:olive')
            plt.plot(*self.slit_wcs.world_to_pixel(bottom), 'o', color='tab:olive')
            plt.plot(*self.slit_wcs.world_to_pixel(top), '^', color='tab:olive')
        plt.show()

    def get_flux(self, slitpos, high=None, low=None, width=None, init_pos=None,
                 allow_reproject=True):
        """

        Parameters
        ----------
        slitpos : SlitParams or astropy.coordinates.SkyCoord
        high : astropy.units.Quantity or None, optional, default is None
            distance between reference coordinates of the given slit
            and its top edge
        low : astropy.units.Quantity or None, optional, default is None
            distance between reference coordinates of the given slit
            and its bottom edge
        width : astropy.units.Quantity or None, optional, default is None
            width of the slit
        init_pos : np.ndarray
            position relative to slit reference coordinate (crpix) on which to
            interpolate the image flux
        allow_reproject : bool, optional, default is True
            if True, reproject if need to (when rotation angle is changed or
            given position is outside the rotated image)

        Returns
        -------

        """
        if type(slitpos) is SlitParams:
            shift = self.slit.refcoord.separation(slitpos.refcoord)
            if (np.abs(slitpos.pa - self.slit.pa) > 0.01
                    or shift.to(u.arcsec).value > 92) and allow_reproject:
                self.rotate_image(slitpos)

            high = slitpos.pos.max() * u.arcsec
            low = slitpos.pos.min() * u.arcsec
            halfw = slitpos.width * u.arcsec / 2.
            x, y = self.slit_wcs.world_to_pixel(slitpos.refcoord)
        else:
            halfw = width / 2.0
            x, y = self.slit_wcs.world_to_pixel(slitpos)

        s_wcs = self.slit_wcs.wcs
        cd_y = s_wcs.cdelt[1] * s_wcs.cunit[1]
        low_pix = (low.to(u.deg) / cd_y.to(u.deg)).value
        high_pix = (high.to(u.deg) / cd_y.to(u.deg)).value
        y_min = int(y + low_pix)
        y_max = int(y + high_pix + 1)
        pos = np.arange(y_max - y_min) + low_pix
        pos = pos * cd_y.to(u.arcsec).value

        cd_x = s_wcs.cdelt[0] * s_wcs.cunit[0]
        hw_x = (halfw.to(u.deg) / cd_x.to(u.deg)).value
        x_min = int(x - hw_x)
        x_max = int(x + hw_x + 1)
        flux = np.sum(self.image_rotated[y_min:y_max, x_min:x_max], axis=1)

        if init_pos is not None:
            flux = np.interp(init_pos, pos, flux)
            pos = init_pos
        else:
            flux = np.interp(slitpos.pos, pos, flux)
            pos = slitpos.pos[:]

        return pos, flux

    @Slot()
    def change_bias(self, bias):
        self.image_rotated += self.bias
        step = (np.max(self.image_hdu.data) - np.min(self.image_hdu.data))/100.
        self.bias = bias * step
        self.image_rotated -= self.bias


class PlotImage:
    """Photometry of an object, its plot and slit shown on image"""

    def __init__(self, figure, frame):
        self.wcs = WCS(frame.header)
        self.figure = figure
        self.axes_obj = figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.image = frame.data
        self.norm_im = simple_norm(self.image, 'linear', percent=99.3)
        self.plot_image()

    def plot_image(self):
        self.axes_obj.clear()
        self.axes_obj = self.figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.axes_obj.imshow(self.image, cmap='bone', norm=self.norm_im)

    def plot_slit(self, slit):
        if self.axes_obj.patches:
            list(self.axes_obj.patches)[1].remove()
            list(self.axes_obj.patches)[0].remove()
            # self.figure.canvas.draw()

        h_up = slit.pos.max() * u.arcsec
        h_down = -slit.pos.min() * u.arcsec
        s = slitPolygon([0 * u.deg, 0 * u.deg], slit.width * u.arcsec, h_up, h_down,
                        theta=0 * u.deg,
                        edgecolor='tab:olive', facecolor='none', lw=0.5,
                        transform=self.axes_obj.get_transform(slit.frame))
        self.axes_obj.add_patch(s)
        t = rotatedTriangle([0 * u.deg, h_up], 5 * slit.width * u.arcsec,
                            theta=0 * u.deg, edgecolor='tab:olive',
                            facecolor='tab:olive', lw=0.5,
                            transform=self.axes_obj.get_transform(slit.frame))
        self.axes_obj.add_patch(t)
        self.figure.canvas.draw()


class PlotSpec:
    def __init__(self, figure, frame, borders=None):
        self.figure = figure
        self.axes_obj = figure.subplots()
        if borders is not None:
            self.x_range = slice(borders[0], borders[1])
            self.y_range = slice(borders[2], borders[3])
        else:
            self.x_range = slice(None, None)
            self.y_range = slice(None, None)
        self.flux = np.sum(frame.data[self.y_range, self.x_range], axis=1)
        self.norm_flux = norm_vector(self.flux)
        self.spectrum = frame.data
        self.hdr = frame.header
        self.pos = self.coords_from_header(axes=2)[self.y_range]
        self.spec_flux_line = None
        self.plot_spectrum_flux()

    def coords_from_header(self, axes=1):
        naxis = self.hdr['NAXIS' + str(axes)]

        if ('CDELT' + str(axes)) in self.hdr:
            cdelt = self.hdr['CDELT' + str(axes)]
        else:
            cdelt = 0.3
        if ('CRPIX' + str(axes)) in self.hdr:
            crpix = self.hdr['CRPIX' + str(axes)]
        else:
            crpix = naxis / 2
        if ('CRVAL' + str(axes)) in self.hdr:
            crval = self.hdr['CRVAL' + str(axes)]
        else:
            crval = 0
        # crval = self.hdr['CRVAL' + str(axes)]
        pix = np.arange(int(naxis)) + 1
        coords = float(cdelt) * (pix - int(crpix)) + float(crval)
        return coords

    def update_scale(self, scale):
        if self.hdr['CDELT2'] == scale:
            return

        self.hdr['CDELT2'] = scale
        self.pos = self.coords_from_header(axes=2)[self.y_range]
        if self.spec_flux_line is not None:
            self.spec_flux_line.set_xdata(self.pos)
            self.axes_obj.relim()
            self.axes_obj.autoscale_view()
            self.figure.canvas.draw()
            self.figure.canvas.flush_events()

    def plot_spectrum_flux(self):
        self.spec_flux_line, = self.axes_obj.plot(self.pos, self.norm_flux, 'b',
                                                  label='spectrum')

    def plot_image_flux(self, pos, flux):
        if len(self.axes_obj.lines) > 1:
            self.axes_obj.lines[1].remove()
        self.axes_obj.plot(pos, norm_vector(flux), 'r', label='image')
        self.axes_obj.legend()
        self.figure.canvas.draw()


class ParameterField(QWidget):

    def __init__(self, spinbox='double', label='',
                 minval: (int, float, Angle, u.Quantity) = 0,
                 maxval: (int, float, Angle, u.Quantity) = 100,
                 dec=2):
        super().__init__()
        self.type = spinbox
        if spinbox == 'double':
            self.spinbox = QDoubleSpinBox()
            self.spinbox.setDecimals(dec)
        elif spinbox == 'int':
            self.spinbox = QSpinBox()
        elif spinbox == 'ra' or spinbox == 'dec':
            self.spinbox = radecSpinBox(radec=spinbox)
        else:
            raise ValueError("spinbox argument of the ParameterField class "
                             + "should be 'int', 'double', 'ra' or 'dec'")

        self.spinbox.setRange(minval, maxval)
        self.spinbox.setKeyboardTracking(False)

        self.label = QLabel(label)
        self.checkbox = QCheckBox('fix')
        self.checkbox.setToolTip("don't change while finding optimal parameters")

        layout = QHBoxLayout()
        layout.addWidget(self.label, 1)
        layout.addWidget(self.spinbox, 6)
        layout.addWidget(self.checkbox, 0)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

    @Slot()
    def setValue(self, val):
        self.spinbox.setValue(val)

    def value(self):
        return self.spinbox.value()

    def getAngle(self):
        if self.type in ['double', 'int']:
            raise AttributeError
        elif self.type in ['ra', 'dec']:
            return self.spinbox.getAngle()

    def getText(self):
        if self.type in ['double', 'int']:
            return str(self.spinbox.value())
        elif self.type in ['ra', 'dec']:
            return self.spinbox.getText()


class PaField(ParameterField):
    def __init__(self, spinbox='double', label='',
                 minval: (int, float, Angle, u.Quantity) = 0,
                 maxval: (int, float, Angle, u.Quantity) = 100,
                 dec=2):
        super().__init__(spinbox, label, minval, maxval, dec)
        self.inverse_shortcut = QShortcut(QKeySequence('Shift+Up'), self)

        self.inverse_shortcut.activated.connect(self.rot_180)

    @Slot()
    def rot_180(self):
        if self.spinbox.hasFocus():
            val = self.spinbox.value()
            val = (val + 180.) % 360.
            self.spinbox.setValue(val)


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    '''
    finished = Signal()  # QtCore.Signal
    error = Signal(tuple)
    result = Signal(object)


class Worker(QRunnable):
    # result = Signal(int)
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @Slot()  # QtCore.Slot
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        res = self.fn(*self.args, **self.kwargs)
        self.signals.result.emit(res)


class PlotWidget(QWidget):
    """main window"""

    def __init__(self, parent=None, spec=None, obj=None, refcenter=None, pa=0.,
                 scale=0.):
        super().__init__(parent)
        self.myStatus = QStatusBar()

        # fits.PrimaryHDU: photometry
        self.image_frame = None
        # fits.PrimaryHdu: spectrum
        self.spec_frame = None
        # PlotImage
        self.image_plot = None
        # PlotSpec
        self.spec_plot = None
        self.interp_params = InterpParams()
        self.slit = SlitParams()
        self.init_slit_frame = self.slit.frame
        self.border_w = SetBorders(self)
        self.bias_w = SetBias(self)
        self.borders = None

        # create widgets
        self.spec_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.image_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.toolbar_spec = NavigationToolbar2QT(self.spec_fig, self)
        self.toolbar_image = NavigationToolbar2QT(self.image_fig, self)

        self.image_field = OpenFile(text='image', mode='o')
        self.spec_field = OpenFile(text='spectrum', mode='o')

        self.PA_input = PaField('double', label='PA, deg',
                                       minval=-360, maxval=360, dec=1)

        self.scale_input = ParameterField('double', label='Scale "/pix',
                                          dec=4)
        # noinspection PyUnresolvedReferences
        self.scale_input.spinbox.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)

        self.ra_input = ParameterField('ra', label='RA, hourangle',
                                       minval=0, maxval=24)
        self.dec_input = ParameterField('dec', label='DEC, degrees',
                                        minval=-90, maxval=90)

        self.x_input = ParameterField('double', label='x, px',
                                      minval=-1e+6, maxval=1e+6)
        self.y_input = ParameterField('double', label='y, px',
                                      minval=-1e+6, maxval=1e+6)
        self.along_input = ParameterField('double',
                                          label='shift along slit, "',
                                          minval=-1e+6, maxval=1e+6)
        self.perp_input = ParameterField('double',
                                         label='shift normal to slit, "',
                                         minval=-1e+6, maxval=1e+6)

        self.inputs_list = [self.PA_input.spinbox, self.scale_input.spinbox,
                            self.ra_input.spinbox, self.dec_input.spinbox,
                            self.x_input.spinbox, self.y_input.spinbox,
                            self.along_input.spinbox, self.perp_input.spinbox]

        self.redraw_button = QPushButton(text='Reopen files')
        self.saveres_button = QPushButton(text='Save Results')
        self.calculate_button = QPushButton(text='Find optimal parameters')
        self.border_button = QPushButton(text='Set borders')
        self.bias_button = QPushButton(text='Set bias')

        self.interp_shortcut = QShortcut(QKeySequence('Alt+I'), self)
        self.progress_box = QMessageBox(QMessageBox.NoIcon,
                                        "Optimal parameters",
                                        "Calculation in progress",
                                        buttons=QMessageBox.NoButton)
        self.progress_box.setStandardButtons(QMessageBox.StandardButton.NoButton)
        # self.progress_box.setText("Searching for optimal parameters")
        self.threadpool = QThreadPool()

        self.__setLayout()

        self.PA_input.setValue(pa)
        self.scale_input.setValue(scale)
        # if filename of an image is set in the terminal
        if obj is not None:
            self.image_field.fill_string(obj)
            self.image_path_changed()
        # if filename of a spectrum is set in the terminal
        if spec is not None:
            self.spec_field.fill_string(spec)
            self.spec_path_changed()

        if refcenter is not None:
            self.ra_input.setValue(refcenter[0])
            self.dec_input.setValue(refcenter[1])

        # self.update_values('eq')

        self.redraw_button.clicked.connect(self.redraw)
        self.calculate_button.clicked.connect(self.find_optimal_parameters)
        self.border_button.clicked.connect(self.set_borders)
        self.bias_button.clicked.connect(self.set_bias)
        self.ra_input.spinbox.valueChanged.connect(lambda: self.update_plots('eq'))
        self.dec_input.spinbox.valueChanged.connect(lambda: self.update_plots('eq'))
        self.x_input.spinbox.valueChanged.connect(lambda: self.update_plots('im'))
        self.y_input.spinbox.valueChanged.connect(lambda: self.update_plots('im'))
        self.along_input.spinbox.valueChanged.connect(lambda: self.update_plots('slit'))
        self.perp_input.spinbox.valueChanged.connect(lambda: self.update_plots('slit'))
        self.PA_input.spinbox.valueChanged.connect(lambda: self.update_plots('pa'))
        self.scale_input.spinbox.valueChanged.connect(lambda: self.update_plots('scale'))
        self.border_w.accepted.connect(self.spec_borders_changed)
        self.saveres_button.clicked.connect(self.save_results)
        self.interp_shortcut.activated.connect(self.plot_rot_image)
        self.bias_w.image_bias.valueChanged.connect(self.set_im_bias)

    def __setLayout(self):
        # Layout
        r_button_layout = QHBoxLayout()
        r_button_layout.addWidget(self.redraw_button)
        r_button_layout.addWidget(self.saveres_button)

        l_button_layout = QHBoxLayout()
        l_button_layout.addWidget(self.calculate_button)
        l_button_layout.addWidget(self.border_button)
        l_button_layout.addWidget(self.bias_button)

        left_layout = QFormLayout()
        left_layout.addRow(self.spec_field)
        left_layout.addRow(self.ra_input)
        left_layout.addRow(self.x_input)
        left_layout.addRow(self.perp_input)
        left_layout.addRow(self.PA_input)
        left_layout.addRow(l_button_layout)
        right_layout = QFormLayout()
        right_layout.addRow(self.image_field)
        right_layout.addRow(self.dec_input)
        right_layout.addRow(self.y_input)
        right_layout.addRow(self.along_input)
        right_layout.addRow(self.scale_input)
        right_layout.addRow(r_button_layout)

        glayout = QGridLayout()
        glayout.addWidget(self.toolbar_spec, 0, 0)
        glayout.addWidget(self.toolbar_image, 0, 1)
        glayout.addWidget(self.spec_fig, 1, 0)
        glayout.addWidget(self.image_fig, 1, 1)
        glayout.addLayout(left_layout, 2, 0)
        glayout.addLayout(right_layout, 2, 1)
        # glayout.addItem(self.myStatus, 3, 0, columnSpan=2)
        glayout.setRowStretch(0, 0)
        glayout.setRowStretch(1, 1)
        glayout.setRowStretch(2, 0)
        self.setLayout(glayout)

    @Slot()
    def plot_rot_image(self):
        self.interp_params.plot_image(self.slit)

    @Slot()
    def image_path_changed(self):
        try:
            print('Opening image ', self.image_field.files)
            self.image_frame = fits.open(self.image_field.files)[0]
        except FileNotFoundError:
            print('IMAGE NOT FOUND')
            self.image_frame = None

    @Slot()
    def spec_path_changed(self):
        try:
            print('Opening spectrum ', self.spec_field.files)
            for i in self.inputs_list:
                i.blockSignals(True)
            self.spec_frame = fits.open(self.spec_field.files)[0]
            self.slit.from_header(self.spec_frame.header)
            self.init_slit_frame = self.slit.refcoord.skyoffset_frame(rotation=self.slit.pa * u.deg)
            self.fill_fileds_from_slit(self.slit)
            self.border_w.add_image(self.spec_frame)
            self.borders = self.border_w.borders
            for i in self.inputs_list:
                i.blockSignals(False)
        except FileNotFoundError:
            print('SPECTRUM NOT FOUND')
            self.spec_frame = None
            for i in self.inputs_list:
                i.blockSignals(False)

    @Slot()
    def spec_borders_changed(self):
        self.borders = self.border_w.borders
        self.slit.add_limits(self.borders)
        try:
            self.image_plot.plot_slit(self.slit)
        except AttributeError:
            print('No object presented')

        if self.spec_frame:
            self.spec_fig.figure.clear()
            self.spec_plot = PlotSpec(self.spec_fig.figure, self.spec_frame,
                                      self.borders)

        if self.image_frame and self.spec_frame:
            self.interp_params.update(self.image_frame, self.slit)
            pos, flux = self.interp_params.get_flux(self.slit,
                                                    init_pos=self.spec_plot.pos)
            self.spec_plot.plot_image_flux(pos, flux)

    @Slot()
    def redraw(self):
        """Redraw all plots. This function runs only when redraw button is clicked."""
        self.spec_path_changed()
        self.image_path_changed()

        # if the string with image frame path is not empty
        if self.image_frame:
            self.image_fig.figure.clear()
            # PlotImage stores image, its wcs and current slit_position
            self.image_plot = PlotImage(self.image_fig.figure, self.image_frame)
            self.image_plot.plot_slit(self.slit)

        if self.spec_frame:
            self.spec_fig.figure.clear()
            self.spec_plot = PlotSpec(self.spec_fig.figure, self.spec_frame)

        if self.image_frame and self.spec_frame:
            self.interp_params.update(self.image_frame, self.slit)
            pos, flux = self.interp_params.get_flux(self.slit,
                                                    init_pos=self.spec_plot.pos)
            self.spec_plot.plot_image_flux(pos, flux)

        self.update_values('eq')
        self.image_fig.draw()
        self.spec_fig.draw()

    def update_values(self, coord=None):
        for i in self.inputs_list:
            i.blockSignals(True)

        if coord == 'eq' or coord == 'pa':
            self.im_from_eq()
            self.slit_from_eq()
        elif coord == 'im':
            self.eq_from_im()
            self.slit_from_eq()
        elif coord == 'slit':
            self.eq_from_slit()
            self.im_from_eq()

        self.slit.from_fields(self.ra_input.getAngle(), self.dec_input.getAngle(),
                              self.PA_input.value(), self.scale_input.value())

        if coord == 'scale' and self.spec_plot is not None:
            self.spec_plot.update_scale(self.scale_input.value())

        for i in self.inputs_list:
            i.blockSignals(False)

    def im_from_eq(self):
        if self.image_plot is not None:
            eq = SkyCoord(self.ra_input.getAngle(), self.dec_input.getAngle(),
                          unit=(u.hourangle, u.deg))
            im = self.image_plot.wcs.world_to_pixel(eq)
            self.x_input.setValue(im[0])
            self.y_input.setValue(im[1])

    def slit_from_eq(self):
        if self.slit is not None:
            eq = SkyCoord(self.ra_input.getAngle(), self.dec_input.getAngle(),
                          unit=(u.hourangle, u.deg))
            slit = eq.transform_to(self.init_slit_frame)
            self.perp_input.setValue(slit.lon.to(u.arcsec).value)
            self.along_input.setValue(slit.lat.to(u.arcsec).value)

    def eq_from_slit(self, coord=None):
        if coord is not None:
            if self.slit is None:
                raise AttributeError('Slit is not set!')
            slit = SkyCoord(coord[0] * u.arcsec,
                            coord[1] * u.arcsec,
                            frame=self.init_slit_frame)
            eq = slit.transform_to(ICRS)
            return eq.ra.to(u.hourangle).value, eq.dec.to(u.deg).value
        if self.slit is not None:
            slit = SkyCoord(self.perp_input.value() * u.arcsec,
                            self.along_input.value() * u.arcsec,
                            frame=self.init_slit_frame)
            eq = slit.transform_to(ICRS)
            self.dec_input.setValue(eq.dec)
            self.ra_input.setValue(eq.ra)

    def eq_from_im(self, coord=None):
        if coord is not None:
            if self.image_frame is None:
                raise AttributeError('Image is not set!')
            eq = self.image_plot.wcs.pixel_to_world(*coord)
            return eq.ra.to(u.hourangle).value, eq.dec.to(u.deg).value
        if self.image_frame is not None:
            im = [self.x_input.value(), self.y_input.value()]
            eq = self.image_plot.wcs.pixel_to_world(*im)
            self.dec_input.setValue(eq.dec)
            self.ra_input.setValue(eq.ra)

    def fill_fileds_from_slit(self, slit: SlitParams):
        for i in self.inputs_list:
            i.blockSignals(True)
        self.PA_input.setValue(slit.pa)
        self.scale_input.setValue(slit.scale)
        self.dec_input.setValue(slit.refcoord.dec)
        self.ra_input.setValue(slit.refcoord.ra)
        self.update_values('eq')
        # self.update_values('scale')
        for i in self.inputs_list:
            i.blockSignals(False)

    @Slot()
    def update_plots(self, coord=None):
        if coord == "pa":
            orig = SkyCoord(self.init_slit_frame.origin)
            rot = self.PA_input.value() * u.deg
            self.init_slit_frame = orig.skyoffset_frame(rotation=rot)
        # Сравниваем значение в поле со значением в интерп (не в слит!)
        need_reproject = (np.abs(self.PA_input.value() - self.interp_params.slit.pa) > 0.01)
        self.update_values(coord)
        if need_reproject:
            self.interp_params.rotate_image(self.slit)

        try:
            self.image_plot.plot_slit(self.slit)
        except AttributeError:
            print('No object presented')

        try:
            pos, flux = self.interp_params.get_flux(self.slit,
                                                    init_pos=self.spec_plot.pos)
            interpflux = np.interp(self.spec_plot.pos, pos, flux)
            print('Qfunc: ', correlation(interpflux, self.spec_plot.flux))
            self.spec_plot.plot_image_flux(pos, flux)
        except AttributeError:
            print("Can't plot flux")

    @Slot()
    def find_optimal_parameters(self):
        self.progress_box.show()
        self.interp_params.rotate_image(self.slit)
        if (self.x_input.checkbox.isChecked()
                or self.y_input.checkbox.isChecked()):
            mode = 'im'
            coords = [self.x_input.value(), self.y_input.value()]
            transform_f = self.eq_from_im

            # TODO: ПОМЕНЯТЬ! ЗАВИСИТ ОТ МАСШТАБА КАРТИНКИ!
            if self.x_input.checkbox.isChecked():
                dra = 0
            else:
                dra = 100

            if self.y_input.checkbox.isChecked():
                ddec = 0
            else:
                ddec = 100
        elif (self.along_input.checkbox.isChecked()
              or self.perp_input.checkbox.isChecked()):
            mode = 'slit'
            coords = [self.along_input.value(), self.perp_input.value()]
            transform_f = self.eq_from_slit

            if self.along_input.checkbox.isChecked():
                dra = 0
            else:
                dra = 30

            if self.perp_input.checkbox.isChecked():
                ddec = 0
            else:
                ddec = 30
        else:
            mode = 'eq'
            coords = [self.slit.refcoord.ra.to(u.hourangle).value,
                      self.slit.refcoord.dec.to(u.deg).value]
            transform_f = lambda x: (x[0], x[1])

            if self.ra_input.checkbox.isChecked():
                dra = 0
            else:
                dra = (30 * u.arcsec).to(u.hourangle).value

            if self.dec_input.checkbox.isChecked():
                ddec = 0
            else:
                ddec = (30 * u.arcsec).to(u.deg).value

        if self.scale_input.checkbox.isChecked():
            dscale = 0
        else:
            dscale = 0.05 * self.slit.scale
        params = [*coords,
                  self.slit.pa,
                  self.slit.scale]

        print(mode)
        print('Params before minimization: ', params)
        bounds = [(params[0] - dra, params[0] + dra),
                  (params[1] - ddec, params[1] + ddec),
                  (params[2], params[2]),
                  (params[3] - dscale, params[3] + dscale)]
        print('Qval before ', self.qfunc_eq(params, self.spec_plot,
                                            self.interp_params, transform_f))
        optargs = (self.spec_plot, self.interp_params, transform_f)
        worker = Worker(self.try_different_minimizers, self.qfunc_eq, params,
                                                    optargs, bounds)
        worker.signals.result.connect(lambda x: self.minimization_finished(x, mode, transform_f))
        # worker.signals.result.connect(print)
        self.threadpool.start(worker)
        # good_params = self.try_different_minimizers(self.qfunc_eq, params,
        #                                             optargs, bounds)

    def minimization_finished(self, good_params, mode, transform_f):
        print(good_params)
        # print(good_params.x)
        print(self.qfunc_eq(good_params.x, self.spec_plot, self.interp_params,
                            transform_f))

        initcoord = SkyCoord(self.slit.refcoord.ra.to(u.hourangle),
                             self.slit.refcoord.dec.to(u.deg))

        for i in self.inputs_list:
            i.blockSignals(True)
        if mode == 'eq':
            self.ra_input.setValue(good_params.x[0] * u.hourangle)
            self.dec_input.setValue(good_params.x[1] * u.deg)
        elif mode == 'im':
            self.x_input.setValue(good_params.x[0])
            self.y_input.setValue(good_params.x[1])
        elif mode == 'slit':
            self.along_input.setValue(good_params.x[0])
            self.perp_input.setValue(good_params.x[1])

        self.PA_input.setValue(good_params.x[2])
        self.scale_input.setValue(good_params.x[3])
        self.update_values('scale')
        for i in self.inputs_list:
            i.blockSignals(False)
        self.calculate_button.setText('Find optimal parameters')
        # self.slit.from_fields(good_params.x[0] * u.hourangle,
        #                       good_params.x[1] * u.deg,
        #                       good_params.x[2], good_params.x[3])
        # self.fill_fileds_from_slit(self.slit)
        # self.update_values('scale')
        self.update_plots(mode)
        rescoord = SkyCoord(self.slit.refcoord.ra.to(u.hourangle),
                            self.slit.refcoord.dec.to(u.deg))
        print('shift is ', initcoord.separation(rescoord))
        self.progress_box.done(0)

    @staticmethod
    def try_different_minimizers(func, params, optargs, bounds, verbose=True):
        """This function use different  methods of scipy.optimize.minimize
        untill they all will consider result as minimum.

        Parameters
        ----------
        func : function to minimize
            will be passed to minimizers
        params : list of float
            initial guess of variable parameters
        optargs : tuple
            additional parameters passed to func
        bounds : list of tuples
            every tuple is the lower and the upper bounds for every parameter
        verbose : bool, optional, default: True
            set to False if you don't want to print information about
            minimization process

        Returns
        -------
        res : result of optimize.minimize

        """
        methods = ['Nelder-Mead', 'Powell', 'L-BFGS-B', 'TNC', 'SLSQP']
        res = None
        q = func(params, *optargs)
        best_params = params
        good_result = False
        while not good_result:
            q_old = q
            for m in methods:
                if verbose: print(m)
                res = optimize.minimize(func, params,
                                        args=optargs, bounds=bounds,
                                        method=m)
                if res.fun < q:
                    q = res.fun
                    if verbose: print('BEST!!!', m)
                    best_params = res.x
            params = best_params
            # if result hasn't changed after checking every minimizer,
            # then consider it's good
            good_result = (q == q_old)
        return res

    @Slot()
    def set_borders(self):
        self.border_w.show()

    @Slot()
    def set_bias(self):
        self.bias_w.show()

    @Slot()
    def set_im_bias(self, bias):
        self.interp_params.change_bias(bias)
        self.update_plots()

    @staticmethod
    def qfunc_eq(params, plotspec, interpparams, transform_f):
        """

        Parameters
        ----------
        params : list of float
            list of parameters to be minimized
            ra - right accession (hourangle)
            dec - declination (degrees)
            pa - position algle (degrees)
            scale - size of pixel ("/pix)
        plotspec : PlotSpec
            object of PlotSpec class, it has flux (summ of a spectrum across
            wavelengths) and position of that flux
        interpparams : InterpParams
            object of InterpParams flux,
            it has photometry roteted to the slit PA
        transform_f : function

        Returns
        -------
        q : float
            this parameter shows how close is the estimation of slit position
            (how good the photometry fits the spectrum flux with given params)

        """
        x, y, pa, scale = params
        ra, dec = transform_f([x, y])

        loc_slit = SlitParams()
        loc_slit.from_header(plotspec.hdr)
        lims = plotspec.y_range
        lims = [None, None, lims.start, lims.stop]
        loc_slit.add_limits(lims)
        loc_slit.set_pos(plotspec.pos)
        loc_slit.from_fields(ra, dec, pa, scale)
        pos, flux = interpparams.get_flux(loc_slit, allow_reproject=True)
        q = correlation(flux, plotspec.flux)
        return q

    @Slot()
    def save_results(self):
        files_path = QFileDialog.getSaveFileName(self, "Save result",
                                                 self.spec_field.dir,
                                                 "All (*)")[0]
        hdul = fits.open(self.spec_field.files)
        hdul[0].header['CDELT2'] = self.scale_input.value()
        hdul[0].header['RA'] = self.ra_input.getText()
        hdul[0].header['DEC'] = self.dec_input.getText()
        if 'EPOCH' in hdul[0].header:
            hdul[0].header['EPOCH'] = 2000.0
        hdul[0].header['POSANG'] = self.PA_input.value()
        if 'CRPIX2' not in hdul[0].header:
            hdul[0].header['CRPIX2'] = self.slit.crpix_init
        hdul.writeto(files_path, overwrite=True)

        print(files_path)
        pass


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--spectrum', nargs='?', default=None,
                        help='''long-slit spectrum''')
    parser.add_argument('-i', '--image', nargs='?', default=None,
                        help='''photometry''')
    pargs = parser.parse_args(sys.argv[1:])

    app = QApplication(sys.argv)
    widg = PlotWidget(None, pargs.spectrum, pargs.image)
    widg.show()
    sys.exit(app.exec())
