from astropy.coordinates import Angle
import astropy.units as u
from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QLineEdit,
    QAbstractSpinBox,
)


class radecSpinBox(QAbstractSpinBox):
    valueChanged = Signal(Angle)

    def __init__(self, parent=None, radec='dec', value=0):
        super().__init__(parent)

        self.line = QLineEdit()
        self.setLineEdit(self.line)

        if radec == 'dec':
            self.unit = u.deg
            self.step = Angle('0:0:1', unit=u.deg)
        elif radec == 'ra':
            self.unit = u.hourangle
            self.step = Angle('0:0:0.1', unit=u.hourangle)

        self.setAccelerated(True)
        self.angle = Angle(value, unit=self.unit)

        self.editingFinished.connect(self.valueFromText)

        self.line.setText(self.textFromValue(self.angle.value))

    def textFromValue(self, val):
        return self.angle.to_string(unit=self.unit, sep=':')

    @Slot()
    def valueFromText(self):
        text = self.text()
        self.angle = Angle(text, unit=self.unit)
        self.line.setText(self.textFromValue(self.angle.value))
        self.valueChanged.emit(self.angle)
        return self.angle.value

    def stepEnabled(self):
        ret = QAbstractSpinBox.StepNone
        ret |= QAbstractSpinBox.StepUpEnabled
        ret |= QAbstractSpinBox.StepDownEnabled
        return ret

    def stepBy(self, steps):
        self.angle += steps * self.step
        self.line.setText(self.textFromValue(self.angle.value))
        self.valueChanged.emit(self.angle)

    def getAngle(self):
        return self.angle

    def setValue(self, value):
        self.angle = Angle(value, unit=self.unit)
        self.line.setText(self.textFromValue(self.angle.value))
        self.valueChanged.emit(self.angle)
