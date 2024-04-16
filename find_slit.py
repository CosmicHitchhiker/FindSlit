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
from astropy.coordinates import SkyCoord, FK5
import reproject
import time
# from itertools import zip_longest, chain

# from matplotlib import pyplot as plt
from matplotlib.figure import Figure
# noinspection PyUnresolvedReferences
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PySide6.QtCore import Slot
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QDoubleSpinBox,
    QHBoxLayout,
    QFormLayout,
    QGridLayout,
    QPushButton,
    QStatusBar,
    # QFileDialog,
    # QCheckBox,
)

from OpenFile import OpenFile
from radecSpinBox import radecSpinBox
from slitPolygon import slitPolygon
matplotlib.use('QtAgg')


def norm_vector(vec):
    return vec / np.linalg.norm(vec)


class SlitParams:
    def __init__(self):
        self.refcoord = SkyCoord(0, 0, unit=(u.hourangle, u.deg))
        self.crpix = 0
        self.npix = 0
        self.PA = 0
        self.scale = 1
        self.width = 1
        self.frame = self.refcoord.skyoffset_frame(rotation=self.PA * u.deg)

    def from_header(self, hdr):
        if 'EPOCH' in hdr:
            equinox = "J" + str(hdr['EPOCH'])
            self.refcoord = SkyCoord(hdr['RA'], hdr['DEC'],
                                     unit=(u.hourangle, u.deg),
                                     frame=FK5(equinox=equinox))
        else:
            self.refcoord = SkyCoord(hdr['RA'], hdr['DEC'],
                                     unit=(u.hourangle, u.deg))

        self.scale = hdr['CDELT2']  # arcsec/pix
        self.npix = hdr['NAXIS2']

        if 'CRPIX2' in hdr:
            self.crpix = hdr['CRPIX2']
        else:
            self.crpix = hdr['NAXIS2'] / 2

        if 'POSANG' in hdr:
            # TDS
            self.PA = hdr['POSANG']
        elif 'PARANGLE' in hdr and 'ROTANGLE' in hdr:
            # SCORPIO
            self.PA = hdr['PARANGLE'] - hdr['ROTANGLE'] + 132.5
        else:
            self.PA = 0

        if self.scale < 0:
            self.PA += 180
            self.scale *= -1

        # if PA not in [0,360) (including negative PA)
        self.PA = self.PA - 360 * (self.PA // 360)

        if 'SLIT' in hdr:
            self.width = self.parse_tds_slit(hdr['SLIT'])
        else:
            self.width = 1.0

        self.frame = self.refcoord.skyoffset_frame(rotation=self.PA * u.deg)

    def from_fields(self, RA, DEC, PA, scale):
        self.PA = PA
        self.scale = scale
        self.refcoord = SkyCoord(RA, DEC,
                                 unit=(u.hourangle, u.deg))
        self.frame = self.refcoord.skyoffset_frame(rotation=self.PA * u.deg)

    def parse_tds_slit(self, slitname):
        if slitname == '1asec':
            return 1.0
        elif slitname == '1.5asec':
            return 1.5
        elif slitname == '10asec':
            return 10
        print('Unknown SLIT keyword')
        return 1.0


class InterpParams:
    """Здесь повёрнутое изображение, можно брать его куски"""
    def __init__(self, image_hdu=None, slit=None):
        self.image_hdu = image_hdu
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

    def get_rot_matrix(self, angle):
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
        slit_PA = slit.PA
        slitpos = slit.refcoord.fk5

        if cdelt is None:
            cdelt = (0.5 * u.arcsec).to(u.deg).value
        if shape is None:
            shape = (500, 1000)
        if center is None:
            center = [shape[0] / 2.0, shape[1] / 2.0]

        self.slit_wcs = WCS(naxis=2)
        self.slit_wcs.wcs.cdelt = [cdelt, cdelt]
        self.slit_wcs.wcs.cunit = [u.deg, u.deg]
        self.slit_wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        self.slit_wcs.wcs.pc = self.get_rot_matrix(slit_PA)
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

        self.slit = slit
        self.slit_hdr = w_header
        return self.slit_wcs, w_header

    def rotate_image(self, slit=None):
        t = time.perf_counter()
        if slit is not None:
            self.slit = slit
            self.make_slit_wcs(slit)
        if self.image_hdu is not None:
            self.image_rotated, _ = reproject.reproject_interp(self.image_hdu,
                                                               self.slit_hdr)
        print('reproject time ', time.perf_counter()-t)

    @Slot()
    def plot_image(self, slitpos=None):
        plt.figure()
        plt.subplot(projection=self.slit_wcs)
        plt.imshow(self.image_rotated)

        if slitpos is not None:
            high = (slitpos.npix - slitpos.crpix) * slitpos.scale * u.arcsec
            low = slitpos.crpix * slitpos.scale * u.arcsec
            halfw = slitpos.width * u.arcsec / 2.
            x, y = self.slit_wcs.world_to_pixel(slitpos.refcoord)
            plt.plot(x, y, 'ro')

        plt.show()

    def get_flux(self, slitpos, high=None, low=None, width=None):
        if type(slitpos) is SlitParams:
            high = (slitpos.npix - slitpos.crpix) * slitpos.scale * u.arcsec
            low = slitpos.crpix * slitpos.scale * u.arcsec
            halfw = slitpos.width * u.arcsec / 2.
            x, y = self.slit_wcs.world_to_pixel(slitpos.refcoord)
        else:
            halfw = width / 2.0
            x, y = self.slit_wcs.world_to_pixel(slitpos)

        s_wcs = self.slit_wcs.wcs
        cd_y = s_wcs.cdelt[1] * s_wcs.cunit[1]
        low_pix = (low.to(u.deg) / cd_y.to(u.deg)).value
        high_pix = (high.to(u.deg) / cd_y.to(u.deg)).value
        y_min = int(y - low_pix)
        y_max = int(y + high_pix + 1)
        pos = np.arange(y_max - y_min) - low_pix
        pos = pos * cd_y.to(u.arcsec)

        cd_x = s_wcs.cdelt[0] * s_wcs.cunit[0]
        hw_x = (halfw.to(u.deg) / cd_x.to(u.deg)).value
        x_min = int(x - hw_x)
        x_max = int(x + hw_x + 1)
        flux = np.sum(self.image_rotated[y_min:y_max, x_min:x_max], axis=1)

        return pos, flux
        

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
            list(self.axes_obj.patches)[0].remove()
            # self.figure.canvas.draw()

        h_up = (slit.npix - slit.crpix) * slit.scale * u.arcsec
        h_down = slit.crpix * slit.scale * u.arcsec
        s = slitPolygon([0 * u.deg, 0 * u.deg], slit.width * u.arcsec, h_up, h_down,
                        theta=0 * u.deg,
                        edgecolor='tab:olive', facecolor='none', lw=0.5,
                        transform=self.axes_obj.get_transform(slit.frame))
        self.axes_obj.add_patch(s)
        self.figure.canvas.draw()


class PlotSpec:
    def __init__(self, figure, frame):
        self.figure = figure
        self.axes_obj = figure.subplots()
        x_range = slice(None, None)
        self.flux = np.sum(frame.data[:, x_range], axis=1)
        self.norm_flux = norm_vector(self.flux)
        self.spectrum = frame.data
        self.hdr = frame.header
        self.pos = self.coords_from_header(axes=2)
        self.spec_flux_line = None
        self.plot_spectrum_flux()

    def coords_from_header(self, axes=1):
        naxis = self.hdr['NAXIS'+str(axes)]
        cdelt = self.hdr['CDELT'+str(axes)]
        crpix = self.hdr['CRPIX'+str(axes)]
        crval = self.hdr['CRVAL'+str(axes)]
        pix = np.arange(int(naxis)) + 1
        coords = float(cdelt) * (pix - int(crpix)) + float(crval)
        return coords

    def update_scale(self, scale):
        if self.hdr['CDELT2'] == scale:
            return

        self.hdr['CDELT2'] = scale
        self.pos = self.coords_from_header(axes=2)
        if self.spec_flux_line is not None:
            self.spec_flux_line.set_xdata(self.pos)
            self.axes_obj.relim()
            self.axes_obj.autoscale_view()
            self.figure.canvas.draw()
            self.figure.canvas.flush_events()

    def plot_spectrum_flux(self):
        self.spec_flux_line,  = self.axes_obj.plot(self.pos, self.norm_flux, 'b')

    def plot_image_flux(self, pos, flux):
        if len(self.axes_obj.lines) > 1:
            self.axes_obj.lines[1].remove()
        self.axes_obj.plot(pos, norm_vector(flux), 'r')
        self.figure.canvas.draw()


class PlotWidget(QWidget):
    """main window"""
    def __init__(self, parent=None, spec=None, obj=None, refcenter=None, PA=0.,
                 scale=0.):
        super().__init__(parent)
        self.myStatus = QStatusBar()

        self.image_frame = None
        self.spec_frame = None
        self.image_plot = None
        self.spec_plot = None
        self.interp_params = InterpParams()
        self.slit = SlitParams()

        # create widgets
        self.spec_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.image_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.toolbar_spec = NavigationToolbar2QT(self.spec_fig, self)
        self.toolbar_image = NavigationToolbar2QT(self.image_fig, self)

        self.image_field = OpenFile(text='image', mode='o')
        self.spec_field = OpenFile(text='spectrum', mode='o')

        self.PA_input = QDoubleSpinBox()
        self.PA_input.setKeyboardTracking(False)
        self.PA_input.setMaximum(360.0)
        self.PA_input.setMinimum(-360.0)

        self.scale_input = QDoubleSpinBox()
        self.scale_input.setKeyboardTracking(False)
        # noinspection PyUnresolvedReferences
        self.scale_input.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)

        self.ra_input = radecSpinBox(radec='ra')
        self.dec_input = radecSpinBox(radec='dec')
        self.ra_input.setKeyboardTracking(False)
        self.dec_input.setKeyboardTracking(False)

        self.redraw_button = QPushButton(text='Reopen files')
        self.saveres_button = QPushButton(text='Save Results')
        self.calculate_button = QPushButton(text='Find optimal parameters')

        # Layout
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.redraw_button)
        button_layout.addWidget(self.saveres_button)

        left_layout = QFormLayout()
        left_layout.addRow(self.spec_field)
        left_layout.addRow('RA', self.ra_input)
        left_layout.addRow('PA', self.PA_input)
        left_layout.addRow(self.calculate_button)
        right_layout = QFormLayout()
        right_layout.addRow(self.image_field)
        right_layout.addRow('DEC', self.dec_input)
        right_layout.addRow('Scale "/pix', self.scale_input)
        right_layout.addRow(button_layout)

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

        self.PA_input.setValue(PA)
        self.scale_input.setValue(scale)
        # if filename of an image is set in the terminal
        if obj is not None:
            self.image_field.fill_string(obj)
            self.imagePathChanged()
        # if filename of a spectrum is set in the terminal
        if spec is not None:
            self.spec_field.fill_string(spec)
            self.specPathChanged()

        if refcenter is not None:
            self.ra_input.setValue(refcenter[0])
            self.dec_input.setValue(refcenter[1])
    #
        self.redraw_button.clicked.connect(self.redraw)
        self.calculate_button.clicked.connect(self.plot_rot_image)
        # self.spec_field.changed_path.connect(self.specChanged)
        # self.image_field.changed_path.connect(self.imagePathChanged)
        self.PA_input.valueChanged.connect(self.updatePlots)
        self.ra_input.valueChanged.connect(self.updatePlots)
        self.dec_input.valueChanged.connect(self.updatePlots)
        self.scale_input.valueChanged.connect(self.updatePlots)
    #     self.ra_input.valueChanged.connect(self.galFrameChanged)
    #     self.dec_input.valueChanged.connect(self.galFrameChanged)
    #     self.vel_input.valueChanged.connect(self.kinematicsChanged)
    #     self.i_input.valueChanged.connect(self.kinematicsChanged)
    #     self.dist_input.valueChanged.connect(self.kinematicsChanged)
    #     self.saveres_button.clicked.connect(self.save_rc)
    #     self.dist_checkbox.stateChanged.connect(self.kinematicsChanged)

    @Slot()
    def plot_rot_image(self):
        self.interp_params.plot_image(self.slit)

    @Slot()
    def imagePathChanged(self):
        try:
            print('Opening image ', self.image_field.files)
            self.image_frame = fits.open(self.image_field.files)[0]
        except FileNotFoundError:
            print('IMAGE NOT FOUND')
            self.image_frame = None

    @Slot()
    def specPathChanged(self):
        try:
            print('Opening spectrum ', self.spec_field.files)
            self.PA_input.blockSignals(True)
            self.ra_input.blockSignals(True)
            self.dec_input.blockSignals(True)
            self.scale_input.blockSignals(True)
            self.spec_frame = fits.open(self.spec_field.files)[0]
            self.slit.from_header(self.spec_frame.header)
            self.fillFiledsFromSlit(self.slit)
            self.PA_input.blockSignals(False)
            self.ra_input.blockSignals(False)
            self.dec_input.blockSignals(False)
            self.scale_input.blockSignals(False)
        except FileNotFoundError:
            print('SPECTRUM NOT FOUND')
            self.spec_frame = None
            self.PA_input.blockSignals(False)
            self.ra_input.blockSignals(False)
            self.dec_input.blockSignals(False)
            self.scale_input.blockSignals(False)

    @Slot()
    def redraw(self):
        """Redraw all plots. This function runs only when redraw button is clicked."""
        self.specPathChanged()
        self.imagePathChanged()

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

        self.image_fig.draw()
        self.spec_fig.draw()

    def updateValues(self):
        self.slit.from_fields(self.ra_input.getAngle(), self.dec_input.getAngle(),
                              self.PA_input.value(), self.scale_input.value())
        self.spec_plot.update_scale(self.scale_input.value())

    def fillFiledsFromSlit(self, slit: SlitParams):
        self.PA_input.setValue(slit.PA)
        self.scale_input.setValue(slit.scale)
        self.dec_input.setValue(slit.refcoord.dec)
        self.ra_input.setValue(slit.refcoord.ra)

    @Slot()
    def updatePlots(self):
        need_reproject = (np.abs(self.PA_input.value() - self.interp_params.slit.PA) > 0.01)
        self.updateValues()
        if need_reproject:
            self.interp_params.rotate_image(self.slit)

        try:
            self.image_plot.plot_slit(self.slit)
            # pos, flux = self.interp_params.get_flux(self.slit)
            # self.spec_plot.plot_image_flux(pos, flux)
        except AttributeError:
            print('No object presented')

        pos, flux = self.interp_params.get_flux(self.slit)
        self.spec_plot.plot_image_flux(pos, flux)
        # try:
        #     pos, flux = self.interp_params.get_flux(self.slit)
        #     # self.spec_plot.plot_image_flux(pos, flux)
        # except AttributeError:
        #     print("Can't plot flux")


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
