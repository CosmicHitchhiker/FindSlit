from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QLineEdit,
    QFileDialog,
    QStyle,
)


class OpenFile(QWidget):
    changed_path = Signal(str)

    def __init__(self, parent=None, text=None, tt=None, mode='n'):
        super().__init__(parent)

        self.fits_box = QLineEdit()
        self.fits_box.setToolTip(tt)
        self._open_folder_action = self.fits_box.addAction(
            qApp.style().standardIcon(QStyle.SP_DirOpenIcon),
            QLineEdit.TrailingPosition)
        self.fits_box.setPlaceholderText(text)
        self.files = None
        self.mode = mode
        self.dir = "/home/astrolander/Documents/Work/DATA"

        vlayout = QVBoxLayout(self)
        vlayout.addWidget(self.fits_box)
        vlayout.setContentsMargins(0, 0, 0, 0)
        vlayout.setSpacing(0)

        self._open_folder_action.triggered.connect(self.on_open_folder)
        self.fits_box.editingFinished.connect(self.check_line)

    @Slot()
    def check_line(self):
        self.files = self.fits_box.text().split(',')
        if self.mode != 'n':
            self.files = self.files[0]
            self.dir = "/".join(self.files.split('/')[:-1])
            self.changed_path.emit(self.dir)
        else:
            self.files = [x.strip() for x in self.files]
            self.dir = "/".join(self.files[0].split('/')[:-1])
            self.changed_path.emit(self.dir)

    @Slot()
    def on_open_folder(self):
        regexps = "All (*)"
        if self.mode == 'n':
            files_path = QFileDialog.getOpenFileNames(self, "Fits", self.dir,
                                                      regexps)[0]
            self.dir = "/".join(files_path[0].split('/')[:-1])
            self.changed_path.emit(self.dir)
        elif self.mode == 'o':
            files_path = QFileDialog.getOpenFileName(self, "Fits", self.dir,
                                                     regexps)[0]
            self.dir = "/".join(files_path.split('/')[:-1])
            self.changed_path.emit(self.dir)
        elif self.mode == 'w':
            files_path = QFileDialog.getSaveFileName(self, "Fits", self.dir,
                                                     regexps)[0]
            self.dir = "/".join(files_path.split('/')[:-1])
            self.changed_path.emit(self.dir)
            ext = files_path.split('.')
            if len(ext) < 2:
                files_path = files_path + '.fits'

        if files_path:
            if self.mode == 'n':
                self.files = files_path.copy()
                self.fits_box.setText(', '.join(files_path))
            else:
                self.files = files_path
                self.fits_box.setText(files_path)

    def fill_string(self, string):
        self.fits_box.setText(string)
        self.check_line()

    def return_filenames(self):
        return self.files
