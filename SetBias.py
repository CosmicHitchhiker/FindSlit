from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QWidget,
    QDialog,
    QDoubleSpinBox,
    QSpinBox,
    QFormLayout,
    QPushButton,
    QLabel,
    QButtonGroup,
)


class SetBias(QDialog):
    newBorders = Signal(list)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.spec_bias = QDoubleSpinBox()
        self.image_bias = QDoubleSpinBox()

        layout = QFormLayout()
        layout.addRow("Spectrum bias", self.spec_bias)
        layout.addRow("Image bias", self.image_bias)
        self.setLayout(layout)