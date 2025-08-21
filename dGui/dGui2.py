# import discamb_py

from RunActionPanel import RunActionPanel
from AtomDisplayStyle import AtomDisplayStyle
from CrystalStructurePresenter import CrystalStructurePresenter
from DisplayPresenter import DisplayPresenter
from DisplayStylePresenter import DisplayStylePresenter
from DisplayStyleComboBox import DisplayStyleComboBox
from FormFactorsModelPresenter import FormFactorsModelPresenter
# from FormFactorsPanel2 import FormFactorsPanel2
from FormFactorsPanel3 import FormFactorsPanel3
from HardwareSettingsDialog import HardwareSettingsDialog
from LassoDrawer import LassoDrawer
from MolCanvas3d import MolCanvas3d
from MolView2 import MolView2
from DisplayControls import DisplayControls
from PathSettingsDialog import PathSettingsDialog
from Settings import Settings
from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6.QtGui import QAction, QIcon, QActionGroup, QFont, QColor, QGuiApplication

import vispy
from vispy.scene import SceneCanvas, visuals
from vispy.app import use_app
from vispy.visuals.filters import ShadingFilter
from scipy.spatial.transform import Rotation
from vispy.visuals.transforms.linear import *
from vispy.visuals.transforms.chain import *
from vispy.scene.visuals import InstancedMesh
from vispy import app, scene, use
from vispy.scene import visuals

import os

use(gl='gl+')

IMAGE_SHAPE = (600, 800)  # (height, width)
NUM_LINE_POINTS = 200
LASSO_COLOR = (1, .1, .1)
PEN_RADIUS = 2
MIN_MOVE_UPDATE_THRESHOLD = 5


class MyMainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.settings = Settings()
        self.setWindowTitle('DiSCaMB GUI')
        my_icon = QIcon()
        my_icon.addFile('img/discamb_icon.bmp')
        self.setWindowIcon(my_icon)
        central_widget = QtWidgets.QWidget()
        main_layout = QtWidgets.QGridLayout()
        self.displayControls = DisplayControls()
        main_layout.addWidget(self.displayControls, 0, 0)

        main_layout.setColumnStretch(0, 1)
        self._canvas3d = MolCanvas3d()
        self._canvas3d.keyPressCallback = self.keyPressEvent
        self._canvas3d.keyReleaseCallback = self.keyReleaseEvent
        self.crystalStructurePresenter = CrystalStructurePresenter(self._canvas3d)
        self.displayPresenter = DisplayPresenter(self.displayControls, self.crystalStructurePresenter)
        main_layout.addWidget(self._canvas3d.canvas.native, 0, 1)

        main_layout.setColumnStretch(1, 10)

        self.drawer = LassoDrawer()
        self.drawer.setVisible(False)
        self.drawer.lassoFinished.connect(self.on_lasso_finish)
        main_layout.addWidget(self.drawer, 0, 1)

        self.tab = QtWidgets.QTabWidget()
        self.tab_model = QtWidgets.QWidget()
        self.tab.addTab(self.tab_model, "Model")
        main_layout.addWidget(self.tab, 0, 2)
        self.formFactorsPanel = FormFactorsPanel3(self.tab_model)
        tab_model_layout = QtWidgets.QVBoxLayout(self.tab_model)
        tab_model_layout.addWidget(self.formFactorsPanel)
        self.formFactorsModelPresenter = FormFactorsModelPresenter(self.formFactorsPanel,
                                                                   self.crystalStructurePresenter)

        self.tab_action = QtWidgets.QWidget()
        self.tab.addTab(self.tab_action, "Run")
        self.actionPanel = RunActionPanel(self.formFactorsModelPresenter, self.tab_model)
        tab_action_layout = QtWidgets.QVBoxLayout(self.tab_action)
        tab_action_layout.addWidget(self.actionPanel)

        # self._controls.lasso_button.clicked.connect(self.switchLasso)
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        self.actionOpen = QAction("&Open", self)
        self.actionOpen.triggered.connect(self.onFileOpen)
        # self.actionOpen.setObjectName(u"actionOpen")
        self.actionClose = QAction("&Close", self)
        # self.actionClose.setObjectName(u"actionClose")
        self.actionClose.triggered.connect(self.onFileClose)
        self.actionExit = QAction("&Exit", self)
        # self.actionExit.setObjectName(u"actionExit")
        self.actionExit.triggered.connect(self.onExit)
        self.actionSelectAll = QAction("Select &All", self)
        self.actionSelectAll.triggered.connect(self.crystalStructurePresenter.selectAll)
        self.actionDeselectAll = QAction("&Deselect All", self)
        self.actionDeselectAll.triggered.connect(self.crystalStructurePresenter.deselect_all)
        self.actionInvertSelection = QAction("&Invert selection")
        self.actionInvertSelection.triggered.connect(self.crystalStructurePresenter.invert_selection)
        self.actionSetPaths = QAction("Set Paths")
        self.actionSetPaths.triggered.connect(self.set_paths)
        self.actionSetHardware = QAction("Set Hardware")
        self.actionSetHardware.triggered.connect(self.set_hardware)
        menubar = self.menuBar()
        menubar.setObjectName(u"menubar")
        menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        file_menu = menubar.addMenu("&File")
        file_menu.addAction(self.actionOpen)
        file_menu.addAction(self.actionClose)
        file_menu.addAction(self.actionExit)
        edit_menu = menubar.addMenu("&Edit")
        edit_menu.addAction(self.actionSelectAll)
        edit_menu.addAction(self.actionDeselectAll)
        edit_menu.addAction(self.actionInvertSelection)
        self.settings_menu = edit_menu.addMenu("Settings")
        self.settings_menu.addAction(self.actionSetPaths)
        self.settings_menu.addAction(self.actionSetHardware)

        display_menu = menubar.addMenu("&Display")
        atom_display_style = display_menu.addMenu("Style")
        style_list = AtomDisplayStyle.get_show_text_list()
        display_style_action_group = QActionGroup(self)
        for style in style_list:
            action = atom_display_style.addAction(style)
            action.setCheckable(True)
            if style == AtomDisplayStyle.default().name.replace(" ", "_"):
                action.setChecked(True)
            display_style_action_group.addAction(action)

        self.display_style_presenter = DisplayStylePresenter(self.displayControls, display_style_action_group,
                                                             self.crystalStructurePresenter)
        label_menu = display_menu.addMenu("&Label")
        self.actionLabelSelected = QAction("&Selected", self)
        self.actionLabelSelected.triggered.connect(self.crystalStructurePresenter.showLabels)
        self.crystalStructurePresenter.showLabels()
        self.actionLabelAll = QAction("&All", self)
        self.actionLabelAll.triggered.connect(self.crystalStructurePresenter.show_all_labels)
        self.actionLabelSelected = QAction("&Selected", self)
        self.actionLabelSelected.triggered.connect(self.crystalStructurePresenter.showLabels)
        self.font = QFont("Helvetica", 12, italic=False)
        self.actionLabelFontStyle = QAction("Set &Font", self)
        self.actionLabelFontStyle.triggered.connect(self.set_font)
        self.actionLabelFontColor = QAction("Set &Color", self)
        self.actionLabelFontColor.triggered.connect(self.set_font_color)
        label_menu.addAction(self.actionLabelSelected)
        label_menu.addAction(self.actionLabelAll)
        label_menu.addAction(self.actionLabelFontStyle)
        label_menu.addAction(self.actionLabelFontColor)

        color_menu = display_menu.addMenu("&Color")
        self.actionSetBackground = QAction("&Background", self)
        self.actionSetBackground.triggered.connect(self.on_set_bcg_color)
        color_menu.addAction(self.actionSetBackground)

        # TOOLBAR

        toolbar = QtWidgets.QToolBar("My main toolbar")
        self.addToolBar(toolbar)

        lasso_icon_file = "img/lasso.png"
        #angle_icon_file = "angle.png"
        distance_icon_file = "img/dist.png"
        label_icon_file = "img/L.png"
        if QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark:
            lasso_icon_file = "img/lasso_dark_mode.png"
            #angle_icon_file = "angle_dark_mode.png"
            distance_icon_file = "img/dist_dark_mode.png"
            label_icon_file = "img/L_dark_mode.png"
        button_lasso = QAction(QIcon(lasso_icon_file), "Lasso selection", self)
        button_lasso.setStatusTip("Lasso selection")
        button_lasso.triggered.connect(self.switchLasso)
        button_lasso.setCheckable(True)
        toolbar.addAction(button_lasso)

        button_labels = QAction(QIcon(label_icon_file), "Show labels", self)
        button_labels.setStatusTip("Show labels")
        button_labels.triggered.connect(self.on_set_labels)
        button_labels.setCheckable(True)
        toolbar.addAction(button_labels)

        button_distance = QAction(QIcon(distance_icon_file), "interatomic distance", self)
        button_distance.setStatusTip("interatomic distance")
        button_distance.triggered.connect(self.crystalStructurePresenter.on_measure_distance_changed)
        button_distance.setCheckable(True)
        toolbar.addAction(button_distance)

        #button_angle = QAction(QIcon(angle_icon_file), "angle", self)
        #button_angle.setStatusTip("angle")
        #button_angle.triggered.connect(self.crystalStructurePresenter.on_measure_angle_changed)
        #button_angle.setCheckable(True)
        #toolbar.addAction(button_angle)

        #toolbar.addWidget(QtWidgets.QLabel(" Labels "))
        #har_toolbar = QtWidgets.QToolBar("HAR toolbar")
        #self.addToolBar(har_toolbar)
        #button_har_plus = QAction(QIcon("plus.png"), "Add chemical unit", self)
        #button_har_plus.setStatusTip("Add chemical unit")
        #button_har_plus.triggered.connect(self.switchLasso)
        #har_toolbar.addAction(button_har_plus)

    def set_paths(self):
        dialog = PathSettingsDialog(files="MATTS")
        dialog.exec()

    def set_hardware(self):
        dialog = HardwareSettingsDialog()
        dialog.hardware_settings_updated.connect(self.formFactorsPanel.harSettingsBox.set_hardware_from_settings)
        dialog.exec()

    def set_font(self):
        font = QtWidgets.QFontDialog.getFont(self.font)
        if font[0]:
            self._canvas3d.labelFont = font[1].family()
            self._canvas3d.labelFontSize = font[1].pointSize()
            self._canvas3d.redrawLabels()

    def set_font_color(self):
        rgb = self._canvas3d.labelColor
        currentColor = QColor(int(rgb[0]), int(rgb[1]), int(rgb[2]))
        color = QtWidgets.QColorDialog.getColor(initial=currentColor)
        if not color.isValid():
            return
        self._canvas3d.labelColor = [color.redF(), color.greenF(), color.blueF()]
        self._canvas3d.redrawLabels()

    def keyPressEvent(self, event):
        modifier = event.keyCombination().keyboardModifiers()
        key = event.key()
        handled = False
        if key == QtCore.Qt.Key_Control:
            if modifier == QtCore.Qt.ControlModifier:
                self.crystalStructurePresenter.setExclusiveSelection(False)
                handled = True
        if key == QtCore.Qt.Key_A:
            if modifier == QtCore.Qt.ControlModifier:
                self.crystalStructurePresenter.selectAll()
                handled = True
        if not handled:
            super().keyPressEvent(event)

    def keyReleaseEvent(self, event):
        if event.key() == QtCore.Qt.Key_Control:
            self.crystalStructurePresenter.setExclusiveSelection(True)
        else:
            super().keyReleaseEvent(event)

    def on_set_labels(self, checked):
        if checked:
            self.crystalStructurePresenter.showLabels()
        else:
            self._canvas3d.unsetLabels()

    def resizeEvent(self, event):
        self._canvas3d.onResize()
        QtWidgets.QMainWindow.resizeEvent(self, event)

    def on_lasso_finish(self, points):
        self._canvas3d.lassoSelection(points)

    def switchLasso(self):
        self.drawer.setVisible(not self.drawer.isVisible())
        # self._canvas_wrapper.view.camera.interactive(not self.drawer.isVisible())

    def onExit(self):
        self.close()

    def onFileClose(self):
        self.crystalStructurePresenter.unset_structure()
        self.formFactorsModelPresenter.unsetCrystalStructure()

    def on_set_bcg_color(self):
        color_rgb = self._canvas3d.canvas.bgcolor.rgb * 255
        color = QColor(int(color_rgb[0]), int(color_rgb[1]), int(color_rgb[2]))
        new_color = QtWidgets.QColorDialog.getColor(QColor.rgb(color), self)
        if new_color.isValid():
            self._canvas3d.set_canvas_color([new_color.red(), new_color.green(), new_color.blue()])

    def onFileOpen(self):
        file_name = QtWidgets.QFileDialog.getOpenFileName(
            self,
            self.tr("Open Crystal Structure File"),
            filter=self.tr("Crystal Structure (*.cif *.ins *.res *.txt)"))
        if file_name[0]:
            self.onFileClose()
            self.formFactorsModelPresenter.project_dir = os.path.dirname(file_name[0])
            self.formFactorsModelPresenter.structure_file = file_name[0]
            self.crystalStructurePresenter.setStructure(file_name[0])


if __name__ == "__main__":
    app = use_app("pyside6")
    app.create()
    s = None
    win = MyMainWindow(s)
    win.show()
    app.run()
