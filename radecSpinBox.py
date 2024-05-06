import numpy as np
from astropy.coordinates import Angle
import astropy.units as u
from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QLineEdit,
    QAbstractSpinBox,
)
from PySide6.QtGui import QValidator


class radecSpinBox(QAbstractSpinBox):
    valueChanged = Signal(Angle)

    def __init__(self, parent=None, radec='dec', value=0):
        super().__init__(parent)

        self.enabled_steps = (QAbstractSpinBox.StepNone
                              | QAbstractSpinBox.StepUpEnabled
                              | QAbstractSpinBox.StepDownEnabled)

        self.line = QLineEdit()
        self.setLineEdit(self.line)

        if radec == 'dec':
            self.unit = u.deg
            self.step = Angle('0:0:1', unit=u.deg)
            self.minimum = Angle(-90, unit=u.deg)
            self.maximum = Angle(90, unit=u.deg)
        elif radec == 'ra':
            self.unit = u.hourangle
            self.step = Angle('0:0:0.1', unit=u.hourangle)
            self.minimum = Angle(0, unit=u.hourangle)
            self.maximum = Angle(24, unit=u.hourangle)

        self.setAccelerated(True)
        self.setKeyboardTracking(False)

        self.precision = 1
        self.dec_precision = self.calc_dec_precision()
        self.angle = Angle(0, unit=self.unit)
        self.setAngle(value)

        self.editingFinished.connect(self.valueFromText)

        self.line.setText(self.textFromValue(self.angle.value))
        self.interpretText()

    def textFromValue(self, val=None):
        if val is None:
            return self.angle.to_string(unit=self.unit, sep=':', precision=self.precision)
        elif isinstance(val, (float, int, str)):
            ang = Angle(val, self.unit)
            return ang.to_string(unit=self.unit, sep=':', precision=self.precision)
        elif isinstance(val, Angle):
            return val.to_string(unit=self.unit, sep=':', precision=self.precision)
        else:
            raise TypeError

    @Slot()
    def valueFromText(self):
        if not self.hasAcceptableInput():
            self.line.setText(self.fixup(self.text()))
        text = self.text()
        self.setAngle(text)
        self.checkBoundaries()
        self.line.setText(self.textFromValue(self.angle.value))
        self.valueChanged.emit(self.angle)
        return self.angle.value

    def stepEnabled(self):
        return self.enabled_steps

    def stepBy(self, steps):
        self.angle += steps * self.step
        self.checkBoundaries()
        self.line.setText(self.textFromValue(self.angle.value))
        self.valueChanged.emit(self.angle)

    def checkBoundaries(self):
        self.enabled_steps = (QAbstractSpinBox.StepNone
                              | QAbstractSpinBox.StepUpEnabled
                              | QAbstractSpinBox.StepDownEnabled)
        if self.angle >= self.maximum:
            self.maximum
            self.enabled_steps ^= QAbstractSpinBox.StepUpEnabled
        if self.angle <= self.minimum:
            self.minimum
            self.enabled_steps ^= QAbstractSpinBox.StepDownEnabled
        self.stepEnabled()
        self.update()
        print(self.enabled_steps)

    def getAngle(self):
        val = round(self.angle.value, self.dec_precision)
        return Angle(val, unit=self.unit)

    def setAngle(self, angle):
        text = self.textFromValue(angle)
        self.angle = Angle(text, unit=self.unit)
        self.line.setText(self.textFromValue(self.angle.value))

    def setValue(self, value):
        self.setAngle(value)
        self.valueChanged.emit(self.angle)

    def validate(self, text, pos):
        """This method is called every time the text in radecSpinBox is changed

        :param text: str, text in the spinbpx
        :param pos: int, coursor position
        :return: QValidator.Acceptable if everythig is ok
                 QValidator.Intermediate if need to call fixup
        """
        # print('validating')
        try:
            Angle(text, unit=self.unit)
        except ValueError:
            # return QValidator.Invalid
            print('need to complete')
            return QValidator.Intermediate
        return QValidator.Acceptable

    def fixup(self, text):
        print('fix')
        return '0:0:0'

    def calc_dec_precision(self):
        minval_sec = 10 ** (-self.precision)
        minval_deg = minval_sec / 3600.
        return int(-np.log10(minval_deg))


