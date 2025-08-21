from PySide6 import QtGui

class ColorItem(QtGui.QStandardItem):
    def __init__(self, color):
        super().__init__()
        self.setEditable(False)
        self.setSelectable(False)
        background = QtGui.QBrush()
        background.setColor(color)
        background.setStyle(QtGui.Qt.BrushStyle.SolidPattern)
        self.setBackground(background)
