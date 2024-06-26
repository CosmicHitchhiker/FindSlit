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
        self.angle = Angle(value, unit=self.unit)

        self.editingFinished.connect(self.valueFromText)

        self.line.setText(self.textFromValue(self.angle.value))
        self.interpretText()

    def textFromValue(self, val=None):
        if val is None:
            return self.angle.to_string(unit=self.unit, sep=':')
        elif isinstance(val, (float, int)):
            ang = Angle(val, self.unit)
            return ang.to_string(unit=self.unit, sep=':')
        elif isinstance(val, Angle):
            return val.to_string(unit=self.unit, sep=':')
        else:
            raise TypeError

    @Slot()
    def valueFromText(self):
        if not self.hasAcceptableInput():
            self.line.setText(self.fixup(self.text()))
        text = self.text()
        self.angle = Angle(text, unit=self.unit)
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
            self.angle = self.maximum
            self.enabled_steps ^= QAbstractSpinBox.StepUpEnabled
        if self.angle <= self.minimum:
            self.angle = self.minimum
            self.enabled_steps ^= QAbstractSpinBox.StepDownEnabled
        self.stepEnabled()
        self.update()

    def getAngle(self):
        return self.angle

    def getText(self):
        return self.line.text()

    def value(self, unit=None):
        """Return angle as float. If no unit is specified, it will be in deg
        for 'dec' spinbox and in hourangle for 'ra' spinbox

        Parameters
        ----------
        unit : astropy.units.Unit
            one of the angle units, if specified, the value will be in this unit

        Returns
        -------
        val : float
            the spinbox value in deg (for 'dec'), hourangle (for 'ra') or
            other specified unit

        """
        if unit:
            return self.angle.to(unit).value
        else:
            return self.angle.to(self.unit).value

    def getUnit(self):
        return self.unit

    def setValue(self, value):
        self.angle = Angle(value, unit=self.unit)
        self.checkBoundaries()
        self.line.setText(self.textFromValue(self.angle.value))
        self.valueChanged.emit(self.angle)

    def setRange(self, minval=None, maxval=None):
        if minval is None:
            pass
        elif isinstance(minval, (Angle, u.Quantity, int, float)):
            self.minimum = Angle(minval, unit=self.unit)
        else:
            raise TypeError('minval should be Angle, Quantity, int, double or None!')
        
        if maxval is None:
            pass
        elif isinstance(maxval, (Angle, u.Quantity, int, float)):
            self.maximum = Angle(maxval, unit=self.unit)
        else:
            raise TypeError('maxval should be Angle, Quantity, int, double or None!')

        self.checkBoundaries()

    def validate(self, text, pos):
        # print('validating')
        try:
            self.angle = Angle(text, unit=self.unit)
        except ValueError:
            # return QValidator.Invalid
            print('need to complete')
            return QValidator.Intermediate
        # return QValidator.Invalid
        return QValidator.Acceptable

    def fixup(self, text):
        print('fix')
        return '0:0:0'
