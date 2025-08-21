from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6.QtGui import QColor
from PySide6.QtGui import QFont
from PySide6.QtGui import QPainter
from PySide6.QtGui import QPen

class LassoDrawer(QtWidgets.QWidget):
    lassoFinished = QtCore.Signal(object)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.x=0
        self.y=0
        self.points=[]
        self.points_last = []
        self.min_move_threshold = 5
        self.last_point = [0,0]
        
    def getPoints(self):
        return self.points_last
    
    
    
    def paintEvent(self, e):
        qp = QPainter()
        qp.begin(self)
        self.drawGeometry(qp)
        qp.end()

    def drawGeometry(self, qp):
        qp.setPen(QPen(QtCore.Qt.green, 4, QtCore.Qt.DashLine))
        #qp.drawEllipse(self.x-30, self.y-30, 60, 60)
        if len(self.points)>2:
          qp.drawPolyline(self.points)  
    
    def mousePressEvent(self, event):
        self.x=event.position().x()
        self.y=event.position().y()
        self.last_point[0]=self.x
        self.last_point[1]=self.y
        self.points.append(QtCore.QPoint(self.x, self.y))    
        self.points.append(QtCore.QPoint(self.x, self.y))
        self.repaint()
        
    def mouseMoveEvent(self, event):
        self.x=event.position().x()
        self.y=event.position().y()
        dx = abs(self.x-self.last_point[0])
        dy = abs(self.y-self.last_point[1])
        if dx>self.min_move_threshold or dy>self.min_move_threshold:
            self.points[-1]=QtCore.QPoint(self.x, self.y)
            self.points.append(self.points[0])    
            self.last_point[0]=self.x
            self.last_point[1]=self.y
            self.repaint()
            
    def mouseReleaseEvent(self, event):        
        self.points_last = [x.toTuple() for x in self.points]
        self.points = []
        self.repaint()
        self.lassoFinished.emit(self.points_last)
        
