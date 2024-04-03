#! /usr/bin/env python3

# Copyright (C) 2022 The Qt Company Ltd.
# SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

import sys

import numpy as np
# from scipy.stats import norm
import matplotlib
from astropy.visualization import simple_norm
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5
# from itertools import zip_longest, chain
from matplotlib import colormaps
from matplotlib import pyplot as plt

# from matplotlib import pyplot as plt
from matplotlib.figure import Figure
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


class SlitParams:
    def __init__(self):
        self.refcoord = SkyCoord(0, 0, unit=(u.hourangle, u.deg))
        self.crpix = 0
        self.npix = 0
        self.PA = 0 * u.deg
        self.scale = 1
        self.width = 1 * u.arcsec
        self.frame = self.refcoord.skyoffset_frame(rotation=self.PA)

    def fromHeader(self, hdr):
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

    def fromFields(self, RA, DEC, PA, scale):
        self.PA = PA * u.deg
        self.scale = scale
        self.refcoord = SkyCoord(RA, DEC,
                                 unit=(u.hourangle, u.deg))
        self.frame = self.refcoord.skyoffset_frame(rotation=self.PA)


    def parse_tds_slit(self, slitname):
        if slitname == '1asec':
            return 1.0
        elif slitname == '1.5asec':
            return 1.5
        elif slitname == '10asec':
            return 10
        print('Unknown SLIT keyword')
        return 1.0


class objectImage():
    '''Photometry of an object, its plot and slit shown on image'''
    def __init__(self, figure, frame):
        self.wcs = WCS(frame.header)
        self.figure = figure
        self.axes_obj = figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.image = frame.data
        self.norm_im = simple_norm(self.image, 'linear', percent=99.3)
        self.plot_image()

    def plot_image(self, gal_frame=None):
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


class PlotSpec():
    def __init__(self, figure, frame):
        self.figure = figure
        self.axes_obj = figure.subplots()
        x_range = slice(None, None)
        self.flux = np.sum(frame.data[:, x_range], axis=1)
        self.spectrum = frame.data
        self.hdr = frame.header
        self.pos = self.coords_from_header(axes=2)
        self.plot_spectrum_flux()

    def coords_from_header(self, axes=1):
        naxis = self.hdr['NAXIS'+str(axes)]
        cdelt = self.hdr['CDELT'+str(axes)]
        crpix = self.hdr['CRPIX'+str(axes)]
        crval = self.hdr['CRVAL'+str(axes)]
        pix = np.arange(int(naxis)) + 1
        coords = float(cdelt) * (pix - int(crpix)) + float(crval)
        return coords

    def plot_spectrum_flux(self):
        self.axes_obj.plot(self.pos, self.flux)



class PlotWidget(QWidget):
    '''main window'''
    def __init__(self, parent=None, spec=None, obj=None, refcenter=None, PA=0.,
                 scale=0.):
        super().__init__(parent)
        self.myStatus = QStatusBar()

        self.image_frame = None
        self.spec_frame = None
        self.image_plot = None
        self.spec_plot = None
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
    def imagePathChanged(self):
        try:
            print('Opening image ', self.image_field.files)
            self.image_frame = fits.open(self.image_field.files)[0]
        except FileNotFoundError:
            print('IMAGE NOT FOUND')
            self.image_frame = None
    #
    @Slot()
    def specPathChanged(self):
        try:
            print('Opening spectrum ', self.spec_field.files)
            self.PA_input.blockSignals(True)
            self.ra_input.blockSignals(True)
            self.dec_input.blockSignals(True)
            self.scale_input.blockSignals(True)
            self.spec_frame = fits.open(self.spec_field.files)[0]
            self.slit.fromHeader(self.spec_frame.header)
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

    #
    #     if self.csv_changed:
    #         self.plot_fig.figure.clear()
    #         data = [pd.read_csv(x) for x in self.csv_field.files]
    #         self.csvGraph = csvPlot(data, self.plot_fig.figure)
    #         slits, masks = self.csvGraph.calc_rc(self.gal_frame,
    #                                              self.inclination,
    #                                              self.sys_vel)
    #         if self.galIm is not None:
    #             self.galIm.plot_galaxy(self.gal_frame)
    #             self.galIm.plot_slit(slits, masks)
    #         self.csv_changed = False
    #
        self.image_fig.draw()
        self.spec_fig.draw()
    #     self.plot_fig.draw()
    #
    # @Slot()
    # def save_rc(self):
    #     filenames = self.csv_field.return_filenames()
    #     dataframes = self.csvGraph.return_rc()
    #     regexps = "CSV (*.csv)"
    #     for fname, dat in zip(filenames, dataframes):
    #         fname_temp = '.'.join(fname.split('.')[:-1]) + '_rc.csv'
    #         file_path = QFileDialog.getSaveFileName(self,
    #                                                 "Save rotation curve",
    #                                                 fname_temp,
    #                                                 regexps)[0]
    #         print('Saving ', file_path)
    #         dat[['RA', 'DEC', 'Circular_v', 'Circular_v_err',
    #              'R_pc', 'R_arcsec', 'mask1', 'mask2']].to_csv(file_path)
    #
    def updateValues(self):
        self.slit.fromFields(self.ra_input.getAngle(), self.dec_input.getAngle(),
                             self.PA_input.value(), self.scale_input.value())
    #     self.inclination = self.i_input.value() * u.deg
    #     self.PA = self.PA_input.value() * u.deg
    #     self.gal_center = SkyCoord(self.ra_input.getAngle(),
    #                                self.dec_input.getAngle(),
    #                                frame='icrs')
    #     self.sys_vel = self.vel_input.value()
    #     self.calc_dist()
    #     self.dist = self.dist_input.value()
    #     self.gal_frame = self.gal_center.skyoffset_frame(rotation=self.PA)

    def fillFiledsFromSlit(self, slit: SlitParams):
        self.PA_input.setValue(slit.PA)
        self.scale_input.setValue(slit.scale)
        self.dec_input.setValue(slit.refcoord.dec)
        self.ra_input.setValue(slit.refcoord.ra)

    @Slot()
    def updatePlots(self):
        self.updateValues()
        try:
            self.image_plot.plot_slit(self.slit)
        except AttributeError:
            print('No object presented')



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--spectrum', nargs='?', default=None,
                        help='''long-slit spectrum''')
    parser.add_argument('-i', '--image', nargs='?', default=None,
                        help='''photometry''')
    pargs = parser.parse_args(sys.argv[1:])

    app = QApplication(sys.argv)
    w = PlotWidget(None, pargs.spectrum, pargs.image)
    w.show()
    sys.exit(app.exec())
