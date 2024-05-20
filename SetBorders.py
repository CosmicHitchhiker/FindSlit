from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QWidget,
    QDialog,
    QDoubleSpinBox,
    QSpinBox,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QButtonGroup,
)
from matplotlib.figure import Figure
# noinspection PyUnresolvedReferences
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from astropy.visualization import simple_norm


class WavelengthTransform:
    def __init__(self, header):
        self.header = header
        self.dl = float(header['CDELT1'])
        self.x0 = float(header['CRPIX1'])
        self.l0 = float(header['CRVAL1'])
        self.x_to_wl = lambda x: (x - self.x0) * self.dl + self.l0
        self.wl_to_x = lambda wl: (wl - self.l0) / self.dl + self.x0

        # print(self.dl, self.x0, self.l0)

        self.lmin = self.x_to_wl(1)
        self.lmax = self.x_to_wl(header['NAXIS1'])


# noinspection PyTypeChecker
class LabledSpinBox(QWidget):
    valueChanged = Signal((int,), (float,))

    def __init__(self, text='', minv=0, maxv=100, value=0, val_type='int'):
        super().__init__()
        if val_type == 'int':
            self.box = QSpinBox()
            self.valueChanged = self.valueChanged[int]
        elif val_type == 'float':
            self.box = QDoubleSpinBox()
            self.valueChanged = self.valueChanged[float]
        else:
            raise AttributeError("'val_type' should be either 'int' or 'float'")
        self.box.setRange(minv, maxv)
        self.box.setValue(value)
        self.label = QLabel(text)

        layout = QHBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.box)
        self.setLayout(layout)

        self.box.valueChanged.connect(self.valueChanged.emit)

    def value(self):
        return self.box.value()

    def set_value(self, value):
        self.box.setValue(value)

    def set_min(self, value):
        self.box.setMinimum(value)

    def set_max(self, value):
        self.box.setMaximum(value)

    def set_range(self, minv, maxv):
        self.box.setRange(minv, maxv)


# noinspection PyTypeChecker
class PairLabledSpinBox(QWidget):
    minChanged = Signal((float,), (int,))
    maxChanged = Signal((float,), (int,))

    def __init__(self, minname='', maxname='', val_type='int'):
        super().__init__()
        self.minbox = LabledSpinBox(minname, val_type=val_type)
        self.maxbox = LabledSpinBox(maxname, val_type=val_type)
        if val_type == 'int':
            self.minChanged = self.minChanged[int]
            self.maxChanged = self.maxChanged[int]
        elif val_type == 'float':
            self.minChanged = self.minChanged[float]
            self.maxChanged = self.maxChanged[float]

        layout = QHBoxLayout()
        layout.addWidget(self.minbox)
        layout.addWidget(self.maxbox)
        self.setLayout(layout)

        self.minbox.valueChanged.connect(self.minChanged)
        self.maxbox.valueChanged.connect(self.maxChanged)
        self.minbox.valueChanged.connect(self.maxbox.set_min)
        self.maxbox.valueChanged.connect(self.minbox.set_max)

    def set_range(self, minv, maxv):
        self.minbox.set_range(minv, maxv)
        self.maxbox.set_range(minv, maxv)

    def set_values(self, minv=None, maxv=None):
        if isinstance(minv, list):
            minv, maxv = minv[0], minv[1]
        if minv is not None:
            self.minbox.set_value(minv)
        if maxv is not None:
            self.maxbox.set_value(maxv)

    def get_min(self):
        return self.minbox.value()

    def get_max(self):
        return self.maxbox.value()


class BorderLine:
    def __init__(self, ax, pos, ltype='h'):
        self.ax = ax
        self.fig = ax.get_figure()
        self.ltype = ltype
        if ltype == 'h':
            self.line = ax.axhline(pos, linestyle='--', color='red', lw=1)
        elif ltype == 'v':
            self.line = ax.axvline(pos, linestyle='--', color='red', lw=1)
        else:
            raise AttributeError("'type' should be 'h' or 'v'")

        self.ltype = ltype

    @Slot()
    def change_pos(self, pos):
        if self.ltype == 'h':
            self.line.set_ydata([pos, pos])
        elif self.ltype == 'v':
            self.line.set_xdata([pos, pos])

        self.fig.canvas.draw()


class FilterButton(QPushButton):
    set_filter = Signal(list)
    def __init__(self, text='', parent=None, limits=None):
        super().__init__(text, parent)
        if limits is not None:
            self.limits = limits
        else:
            self.limits = [3500., 8000.]

        self.clicked.connect(lambda: self.set_filter.emit(self.limits))


class SetBorders(QDialog):
    newBorders = Signal(list)

    def __init__(self, parent=None):
        super().__init__(parent)

        # block any actions in the parent window
        # noinspection PyUnresolvedReferences
        self.setModal(QDialog.Accepted)

        self.fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.toolbar = NavigationToolbar2QT(self.fig, self)
        self.ax = self.fig.figure.subplots()
        self.apply_button = QPushButton(text='Apply')
        self.cancel_button = QPushButton(text='Cancel')
        self.img = None
        self.norm_im = None
        self.cut_img = None
        self.wl_transform = None
        self.borders = []
        self.minx_line = BorderLine(self.ax, 0, 'v')
        self.maxx_line = BorderLine(self.ax, 0, 'v')
        self.miny_line = BorderLine(self.ax, 0, 'h')
        self.maxy_line = BorderLine(self.ax, 0, 'h')
        self.counter = 0

        self.x_box = PairLabledSpinBox('x_min', 'x_max')
        self.y_box = PairLabledSpinBox('y_min', 'y_max')
        self.l_box = PairLabledSpinBox('l_min', 'l_max',
                                       val_type='float')
        filters = {"u": [3000, 4000],
                   "g": [4000, 5350],
                   "r": [5500, 6800],
                   "i": [6900, 8200]}

        self.filter_buttons = dict()
        filt_layout = QHBoxLayout()
        for f in filters:
            self.filter_buttons[f] = FilterButton(f, self, filters[f])
            filt_layout.addWidget(self.filter_buttons[f])

        val_layout = QHBoxLayout()
        val_layout.addWidget(self.x_box)
        val_layout.addWidget(self.y_box)
        val_layout.addWidget(self.l_box)

        b_layout = QHBoxLayout()
        b_layout.addWidget(self.cancel_button)
        b_layout.addWidget(self.apply_button)
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.fig)
        layout.addLayout(filt_layout)
        layout.addLayout(val_layout)
        layout.addLayout(b_layout)
        layout.setStretch(1, 1)
        self.setLayout(layout)

        self.x_box.minChanged.connect(self.minx_changed)
        self.x_box.maxChanged.connect(self.maxx_changed)
        self.y_box.minChanged.connect(self.miny_line.change_pos)
        self.y_box.maxChanged.connect(self.maxy_line.change_pos)
        self.l_box.minChanged.connect(self.minl_changed)
        self.l_box.maxChanged.connect(self.maxl_changed)
        self.cancel_button.clicked.connect(self.reject)
        self.apply_button.clicked.connect(self.apply_borders)
        for f in self.filter_buttons.values():
            f.set_filter.connect(self.l_box.set_values)

    def add_image(self, frame):
        self.img = frame.data
        self.norm_im = simple_norm(self.img, 'linear', percent=85.)
        self.ax.imshow(self.img, cmap='bone', norm=self.norm_im, origin='lower')

        self.wl_transform = WavelengthTransform(frame.header)
        minx = 1
        miny = 1
        maxx = frame.header['NAXIS1']
        maxy = frame.header['NAXIS2']
        minl = self.wl_transform.x_to_wl(minx)
        maxl = self.wl_transform.x_to_wl(maxx)

        self.x_box.set_range(minx, maxx)
        self.y_box.set_range(miny, maxy)
        self.l_box.set_range(minl, maxl)

        self.x_box.set_values(minx, maxx)
        self.y_box.set_values(miny, maxy)
        self.l_box.set_values(minl, maxl)

        # self.minx_line = BorderLine(self.ax, minx-1, 'v')
        self.minx_line.change_pos(minx-1)
        self.maxx_line.change_pos(maxx-1)
        self.miny_line.change_pos(miny-1)
        self.maxy_line.change_pos(maxy-1)

        self.apply_borders()

        self.fig.draw()

    def apply_borders(self):
        self.borders = [self.x_box.get_min() - 1,
                        self.x_box.get_max(),
                        self.y_box.get_min() - 1,
                        self.y_box.get_max()]

        super().accept()

    @Slot()
    def minl_changed(self, val):
        minx = self.wl_transform.wl_to_x(val)
        self.x_box.minbox.set_value(round(minx))

    @Slot()
    def maxl_changed(self, val):
        maxx = self.wl_transform.wl_to_x(val)
        self.x_box.maxbox.set_value(round(maxx))

    @Slot()
    def minx_changed(self, val):
        self.minx_line.change_pos(val)
        self.l_box.blockSignals(True)
        minl = self.wl_transform.x_to_wl(val)
        self.l_box.minbox.set_value(minl)
        self.l_box.blockSignals(False)

    @Slot()
    def maxx_changed(self, val):
        self.maxx_line.change_pos(val)
        self.l_box.blockSignals(True)
        maxl = self.wl_transform.x_to_wl(val)
        self.l_box.maxbox.set_value(maxl)
        self.l_box.blockSignals(False)


