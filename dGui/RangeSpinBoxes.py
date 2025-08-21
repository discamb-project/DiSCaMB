from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6 import QtGui
from PySide6.QtCore import Signal

class RangeSpinBoxes(QtWidgets.QHBoxLayout): 
    changed = Signal(float, float)
    
    def __init__(self, low=0, high=1, minVal=-1e+10, maxVal=1e+10, parent=None):
        super().__init__(parent)
        self.lowSpinBox = QtWidgets.QDoubleSpinBox()
        self.highSpinBox = QtWidgets.QDoubleSpinBox()
        self.low = low
        self.high = high
        self.maxVal = maxVal
        self.minVal = minVal
        self.lowSpinBox.setMaximum(self.high)
        self.highSpinBox.setMaximum(self.maxVal)
        self.lowSpinBox.setMinimum(self.minVal)
        self.highSpinBox.setMinimum(self.low)
        
        self.lowSpinBox.setValue(self.low)
        self.highSpinBox.setValue(self.high)
        self.lowSpinBox.setSingleStep(0.1)
        self.lowSpinBox.valueChanged.connect(self.__onMinChanged)
        self.highSpinBox.valueChanged.connect(self.__onMaxChanged)
        self.highSpinBox.setSingleStep(0.1)
        self.addWidget(self.lowSpinBox)
        self.addWidget(self.highSpinBox)
    
    def getValues(self):
        return self.lowSpinBox.value(), self.highSpinBox.value()
    
    def getMinValue(self):
        return self.lowSpinBox.value()
    
    def getMaxValue(self):
        return self.highSpinBox.value()
        
    def show(self, showWidgets):    
        if showWidgets:
            self.lowSpinBox.show()
            self.highSpinBox.show()
        else:
            self.lowSpinBox.hide()
            self.highSpinBox.hide()
    
    def __onMinChanged(self, value):
        self.low = value
        self.changed.emit(self.low, self.high)
        self.highSpinBox.setMinimum(self.low)

    def __onMaxChanged(self, value):
        self.high = value
        self.changed.emit(self.low, self.high)
        self.lowSpinBox.setMaximum(self.high)
