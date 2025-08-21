from CrystalStructurePresenter import CrystalStructurePresenter
from DisplayControls import DisplayControls
from AtomDisplayStyle import AtomDisplayStyle
from PySide6 import QtGui

from enum import IntEnum


class DisplayStylePresenter:
    def __init__(
            self,
            display_controls: DisplayControls,
            display_style_action_group: QtGui.QActionGroup,
            crystal_presenter: CrystalStructurePresenter):

        self.display_controls = display_controls
        self.display_style_action_group = display_style_action_group
        self.crystal_presenter = crystal_presenter
        self.current_style = AtomDisplayStyle.default()
        self.__init_current_style()
        self.display_style_action_group.triggered.connect(self.on_menu_selection_changed)
        self.display_controls.displayStyleComboBox.currentIndexChanged.connect(self.on_display_box_selection_changed)

    def __set_action_group_style(self, style):
        actions = self.display_style_action_group.actions()
        for action in actions:
            action_style = AtomDisplayStyle[action.text().replace(" ", "_")]
            if action_style == style:
                if not action.isChecked():
                    action.setChecked(True)

    def __init_current_style(self):
        self.current_style = AtomDisplayStyle.default()
        self.crystal_presenter.set_atom_style(self.current_style)
        self.__set_display_style_combo_box_style(self.current_style)
        self.__set_action_group_style(self.current_style)

    def __set_display_style_combo_box_style(self, style):
        if self.display_controls.displayStyleComboBox.currentIndex() == int(style):
            return
        self.display_controls.displayStyleComboBox.set_atom_style(style)

    def set_style(self, style):
        if self.current_style == style:
            return
        self.current_style = style
        self.__set_action_group_style(style)
        self.__set_display_style_combo_box_style(style)
        self.crystal_presenter.set_atom_style(style)

    def on_menu_selection_changed(self, action):
        style_menu = AtomDisplayStyle[action.text().replace(" ", "_")]
        self.set_style(style_menu)

    def on_display_box_selection_changed(self, idx):
        if idx == self.current_style:
            return
        self.set_style(AtomDisplayStyle(self.display_controls.displayStyleComboBox.currentIndex()))


