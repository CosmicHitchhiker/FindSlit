from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QWidget,
    QDialog,
    QDoubleSpinBox,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
)
from matplotlib.figure import Figure
# noinspection PyUnresolvedReferences
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from astropy.visualization import simple_norm


class SetBorders(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        # block any actions in the parent window
        self.setModal(QDialog.Accepted)

        self.fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.toolbar = NavigationToolbar2QT(self.fig, self)
        self.ax = self.fig.figure.subplots()
        self.apply_button = QPushButton(text='Apply')
        self.cancel_button = QPushButton(text='Cancel')
        self.img = None
        self.norm_im = None
        self.counter = 0

        h_layout = QHBoxLayout()
        h_layout.addWidget(self.cancel_button)
        h_layout.addWidget(self.apply_button)
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.fig)
        layout.addLayout(h_layout)
        self.setLayout(layout)

        self.cancel_button.clicked.connect(self.increase_counter)

    @Slot()
    def increase_counter(self):
        self.counter += 1

    def add_image(self, img):
        self.img = img
        self.norm_im = simple_norm(img, 'linear', percent=85.)
        self.ax.imshow(img, cmap='bone', norm = self.norm_im, origin='lower')
        self.fig.draw()

