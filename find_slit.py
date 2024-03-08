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
matplotlib.use('QtAgg')


class slitParams:
    def __init__(self):
        self.refcoord = SkyCoord(0, 0, unit=(u.hourangle, u.deg))
        self.crpix = 0
        self.npix = 0
        self.PA = 0 * u.deg
        self.scale = 1
        self.width = 1 * u.arcsec

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

    def fromFields(self, RA, DEC, PA, scale):
        self.PA = PA
        self.scale = scale
        self.refcoord = SkyCoord(RA, DEC,
                                 unit=(u.hourangle, u.deg))

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
        # self.slits = None
        # self.masks = None
        # self.slit_draws = None
        self.plot_image()

    def plot_image(self, gal_frame=None):
        self.axes_obj.clear()
        self.axes_obj = self.figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.axes_obj.imshow(self.image, cmap='bone', norm=self.norm_im)
        # "Стираем" чёрточки по краям картинки
        # self.axes_gal.coords['ra'].set_ticks(color='white')
        # self.axes_gal.coords['dec'].set_ticks(color='white')

        # if self.slits is not None:
        #     self.plot_slit(self.slits, self.masks)

    def plot_slit(self, slits, masks):
        self.slits = slits
        self.masks = masks
        for slit, mask in zip(slits, masks):
            plot_slit_points(self.axes_gal, slit, mask,
                             'icrs')

        for line in self.axes_gal.lines:
            self.axes_gal.lines.remove(line)

        for slit, mask, i in zip(slits, masks, range(0, 20, 2)):
            mask1, mask2 = mask
            if len(mask1[mask1]) > 0:
                self.axes_gal.plot(
                    slit.ra[mask1],
                    slit.dec[mask1],
                    marker='.',
                    linestyle='',
                    transform=self.axes_gal.get_transform('icrs'),
                    color=self.colors[i])
            if len(mask2[mask2]) > 0:
                self.axes_gal.plot(
                    slit.ra[mask2],
                    slit.dec[mask2],
                    marker='.',
                    linestyle='',
                    transform=self.axes_gal.get_transform('icrs'),
                    color=self.colors[i + 1])


class csvPlot():
    def __init__(self, data, figure):
        self.colors = colormaps['tab20'](np.linspace(0, 1, 20))
        # data - list of pd.DataFrame
        self.data = data
        self.slits = []
        for dat in self.data:
            slit_ra = dat['RA']
            slit_dec = dat['DEC']
            self.slits.append(SkyCoord(slit_ra, slit_dec, frame='icrs',
                                       unit=(u.hourangle, u.deg)))
        self.axes_plot = figure.subplots()

    def calc_rc(self, gal_frame, inclination, sys_vel, dist=None):
        if dist is None:
            dist = sys_vel / 70.
        self.axes_plot.clear()
        self.masks = []
        for dat, slit in zip(self.data, self.slits):
            dat = los_to_rc(dat, slit, gal_frame, inclination, sys_vel, dist)
            self.masks.append([dat['mask1'].to_numpy(),
                               dat['mask2'].to_numpy()])
        self.plot_rc()
        return self.slits, self.masks

    def plot_rc(self):
        self.axes_plot.set_ylabel('Circular Velocity, km/s')
        self.axes_plot.set_xlabel('R, parsec')
        for dat, mask, i in zip(self.data, self.masks, range(0, 20, 2)):
            verr = dat['Circular_v_err'].to_numpy()
            mask1, mask2 = mask
            if len(mask1[mask1]) > 0:
                self.axes_plot.errorbar(
                    dat['R_pc'][mask1],
                    dat['Circular_v'][mask1],
                    yerr=verr[mask1],
                    linestyle='',
                    marker='.',
                    color=self.colors[i])
            if len(mask2[mask2]) > 0:
                self.axes_plot.errorbar(
                    dat['R_pc'][mask2],
                    dat['Circular_v'][mask2],
                    yerr=verr[mask2],
                    linestyle='',
                    marker='.',
                    color=self.colors[i + 1])

    def return_rc(self):
        return self.data


class PlotWidget(QWidget):
    '''main window'''
    def __init__(self, parent=None, spec=None, frame=None, refcenter=None, PA=0.,
                 scale=0., velocity=0.):
        super().__init__(parent)
        self.myStatus = QStatusBar()

        self.obj_frame = None
        self.spec_frame = None
        self.slit = slitParams()

        # flags to replot image and flux
        # self.image_changed = False
        # self.spec_changed = False

        # create widgets
        self.spec_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.image_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.toolbar_spec = NavigationToolbar2QT(self.spec_fig, self)
        self.toolbar_image = NavigationToolbar2QT(self.image_fig, self)

        self.image_field = OpenFile(text='image', mode='o')
        # if filename of an image is set in the terminal
        if frame is not None:
            self.image_field.fill_string(frame)
            # need to plot image
            self.image_changed = True

        self.spec_field = OpenFile(text='spectrum', mode='o')
        # if filename of a spectrum is set in the terminal
        if spec is not None:
            self.spec_field.fill_string(spec)
            self.spec_changed = True

        self.PA_input = QDoubleSpinBox()
        self.PA_input.setKeyboardTracking(False)
        self.PA_input.setMaximum(360.0)
        self.PA_input.setMinimum(-360.0)
        self.PA_input.setValue(PA)

        self.scale_input = QDoubleSpinBox()
        self.scale_input.setKeyboardTracking(False)
        self.scale_input.setValue(scale)

        self.ra_input = radecSpinBox(radec='ra')
        self.dec_input = radecSpinBox(radec='dec')
        self.ra_input.setKeyboardTracking(False)
        self.dec_input.setKeyboardTracking(False)
        if refcenter is not None:
            self.ra_input.setValue(refcenter[0])
            self.dec_input.setValue(refcenter[1])

        self.redraw_button = QPushButton(text='Redraw')
        self.saveres_button = QPushButton(text='Save Results')
        self.calculate_button = QPushButton(text='Find optimal parameters')

        # Layout
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.redraw_button)
        button_layout.addWidget(self.saveres_button)

        left_layout = QFormLayout()
        left_layout.addRow(self.spec_field)
        left_layout.addRow('RA', self.ra_input)
        left_layout.addRow('DEC', self.dec_input)
        left_layout.addRow(self.calculate_button)
        right_layout = QFormLayout()
        right_layout.addRow(self.image_field)
        right_layout.addRow('PA', self.PA_input)
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
        # glayout.setRowStretch(3, 0)
        self.setLayout(glayout)
    #
        self.redraw_button.clicked.connect(self.redraw)
        self.spec_field.changed_path.connect(self.specChanged)
        self.image_field.changed_path.connect(self.imagePathChanged)
    #     self.PA_input.valueChanged.connect(self.galFrameChanged)
    #     self.ra_input.valueChanged.connect(self.galFrameChanged)
    #     self.dec_input.valueChanged.connect(self.galFrameChanged)
    #     self.vel_input.valueChanged.connect(self.kinematicsChanged)
    #     self.i_input.valueChanged.connect(self.kinematicsChanged)
    #     self.dist_input.valueChanged.connect(self.kinematicsChanged)
    #     self.saveres_button.clicked.connect(self.save_rc)
    #     self.dist_checkbox.stateChanged.connect(self.kinematicsChanged)
    #
    @Slot()
    def imagePathChanged(self):
        try:
            print('Opening ', self.image_field.files)
            self.obj_frame = fits.open(self.image_field.files)[0]
        except FileNotFoundError:
            print('IMAGE NOT FOUND')
            self.obj_frame = None
    #
    @Slot()
    def specChanged(self):
        try:
            print('Opening ', self.spec_field.files)
            self.spec_frame = fits.open(self.spec_field.files)[0]
            self.slit.fromHeader(self.spec_frame.header)
            self.fillFiledsFromSlit(self.slit)
        except FileNotFoundError:
            print('SPECTRUM NOT FOUND')
            self.spec_frame = None
    #
    # @Slot()
    # def calc_dist(self):
    #     if self.dist_checkbox.isChecked():
    #         self.dist_input.setValue(self.vel_input.value() / 70.)
    #         self.dist_input.setDisabled(True)
    #     else:
    #         self.dist_input.setDisabled(False)
    #
    # @Slot()
    # def galFrameChanged(self):
    #     self.updateValues()
    #     self.galIm.plot_galaxy(self.gal_frame)
    #     slits, masks = self.csvGraph.calc_rc(self.gal_frame, self.inclination,
    #                                          self.sys_vel)
    #     self.galIm.plot_slit(slits, masks)
    #     self.gal_fig.draw()
    #     self.plot_fig.draw()
    #
    # @Slot()
    # def kinematicsChanged(self):
    #     self.updateValues()
    #     self.csvGraph.calc_rc(self.gal_frame, self.inclination, self.sys_vel,
    #                           self.dist)
    #     self.plot_fig.draw()
    #
    @Slot()
    def redraw(self):
        """Redraw all plots. This function runs only when redraw button is clicked."""
        self.updateValues()

        # if the string with image frame path is not empty
        if self.obj_frame:
            self.image_fig.figure.clear()
            # objectImage stores image, its wcs and current slit_position
            self.objIm = objectImage(self.image_fig.figure, self.obj_frame)
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
        pass
    #     self.inclination = self.i_input.value() * u.deg
    #     self.PA = self.PA_input.value() * u.deg
    #     self.gal_center = SkyCoord(self.ra_input.getAngle(),
    #                                self.dec_input.getAngle(),
    #                                frame='icrs')
    #     self.sys_vel = self.vel_input.value()
    #     self.calc_dist()
    #     self.dist = self.dist_input.value()
    #     self.gal_frame = self.gal_center.skyoffset_frame(rotation=self.PA)

    def fillFiledsFromSlit(self, slit: slitParams):
        self.PA_input.setValue(slit.PA)
        self.scale_input.setValue(slit.scale)
        self.dec_input.setValue(slit.refcoord.dec)
        self.ra_input.setValue(slit.refcoord.ra)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--csv', nargs='+', default=None,
                        help='''csv or similar file with positions and
                        velocities''')
    parser.add_argument('-r', '--refcenter', nargs=2, default=None,
                        help='''coordinates of center of galaxy''')
    parser.add_argument('-v', '--velocity', type=float, default=0.0,
                        help='system velocity')
    parser.add_argument('-p', '--PA', type=float, default=0.0,
                        help='galaxy PA')
    parser.add_argument('-i', '--inclination', type=float, default=0.0,
                        help='inclination of galaxy')
    parser.add_argument('-f', '--frame', default=None,
                        help='frame with image')
    pargs = parser.parse_args(sys.argv[1:])

    app = QApplication(sys.argv)
    w = PlotWidget(None, pargs.csv, pargs.frame, pargs.refcenter, pargs.PA,
                   pargs.inclination, pargs.velocity)
    w.show()
    sys.exit(app.exec())
