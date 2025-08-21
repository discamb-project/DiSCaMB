import copy
from PySide6 import QtWidgets
from AtomDisplayStyle import AtomDisplayStyle


class DisplayStyleComboBox(QtWidgets.QComboBox):

    def __init__(self, atom_display_style, parent=None):
        super().__init__(parent)
        self.addItems(AtomDisplayStyle.get_show_text_list())
        self.setCurrentIndex(int(AtomDisplayStyle.default()))

    def set_atom_style(self, atom_display_style):
        if self.currentIndex() == int(atom_display_style):
            return
        else:
            self.setCurrentIndex(int(atom_display_style))




