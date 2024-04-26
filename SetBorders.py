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
)
from matplotlib.figure import Figure
# noinspection PyUnresolvedReferences
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from astropy.visualization import simple_norm


class WavelengthTransform:
    def __init__(self, header):
        self.header = header
        self.dl = header['CDELT1']
        self.x0 = header['CRPIX1']
        self.l0 = header['CRVAL1']
        self.x_to_wl = lambda x: (x - self.x0) * self.dl + self.l0
        self.wl_to_x = lambda wl: (wl - self.l0) / self.dl + self.x0

        self.lmin = self.x_to_wl(1)
        self.lmax = self.x_to_wl(header['NAXIS1'])


class LabledSpinBox(QWidget):
    valueChanged = Signal((int,),(float,))

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

    def set_values(self, minv, maxv):
        self.minbox.set_value(minv)
        self.maxbox.set_value(maxv)


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



class SetBorders(QDialog):
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
        self.wl_transform = None
        self.minx_line = BorderLine(self.ax, 0, 'v')
        self.maxx_line = BorderLine(self.ax, 0, 'v')
        self.miny_line = BorderLine(self.ax, 0, 'h')
        self.maxy_line = BorderLine(self.ax, 0, 'h')
        self.counter = 0

        self.x_box = PairLabledSpinBox('x_min', 'x_max')
        self.y_box = PairLabledSpinBox('y_min', 'y_max')
        self.l_box = PairLabledSpinBox('l_min', 'l_max',
                                       val_type='float')

        # self.min_x = LabledSpinBox('x_min')
        # self.max_x = LabledSpinBox('x_max')
        # self.min_y = LabledSpinBox('y_min')
        # self.max_y = LabledSpinBox('y_max')
        # self.min_l = LabledSpinBox('l_min')
        # self.max_l = LabledSpinBox('l_max')
        # self.spinboxes = [self.min_x, self.max_x, self.min_y,
        #                   self.max_y, self.min_l, self.max_l]
        # [x.box.setDecimals(0) for x in self.spinboxes[:4]]

        val_layout = QHBoxLayout()
        val_layout.addWidget(self.x_box)
        val_layout.addWidget(self.y_box)
        val_layout.addWidget(self.l_box)
        # [val_layout.addWidget(x) for x in self.spinboxes]

        b_layout = QHBoxLayout()
        b_layout.addWidget(self.cancel_button)
        b_layout.addWidget(self.apply_button)
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.fig)
        layout.addLayout(val_layout)
        layout.addLayout(b_layout)
        self.setLayout(layout)

        self.x_box.minChanged.connect(self.minx_line.change_pos)
        self.x_box.maxChanged.connect(self.maxx_line.change_pos)
        self.y_box.minChanged.connect(self.miny_line.change_pos)
        self.y_box.maxChanged.connect(self.maxy_line.change_pos)

    @Slot()
    def increase_counter(self, num=None):
        if num is not None:
            print(num)
        self.counter += 1

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



        # for box, val in zip(self.spinboxes, vals):
        #     box.set_value(val)

        self.fig.draw()
