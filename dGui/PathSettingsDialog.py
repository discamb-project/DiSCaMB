from Settings import Settings
from PySide6 import QtWidgets
from PySide6.QtCore import Signal


class PathSettingsDialogButton(QtWidgets.QPushButton):
    button_clicked = Signal(str)

    def __init__(self, name, text, parent=None):
        super().__init__(text, parent=parent)
        self.name = name
        super().clicked.connect(self.on_clicked)

    def on_clicked(self):
        self.button_clicked.emit(self.name)


class PathSettingsDialog(QtWidgets.QDialog):

    def __init__(self, files):
        super().__init__()
        self.settings = Settings()
        self.files = files
        layout = QtWidgets.QGridLayout()
        all_keys = self.settings.keys()
        self.path_edits = {}
        self.paths = {}
        idx = 0
        for key in all_keys:
            words = key.split('/')
            if words[0] == 'paths':
                value = ""
                if key in files:
                    value = self.settings.get(key)
                    if value:
                        value = value[0]
                else:
                    value = self.settings.get(key)
                self.paths[words[1]] = value
                label = QtWidgets.QLabel(words[1])
                value_str = str(value)
                if value is None:
                    value_str = ""
                file_name_edit = QtWidgets.QLineEdit(value_str)
                file_name_edit.setReadOnly(True)
                self.path_edits[words[1]] = file_name_edit
                if value is not None:
                    #if isinstance(value, tuple):
                    #    file_name_edit.setText(value[0])
                    #else:
                    file_name_edit.setText(value)
                button = PathSettingsDialogButton(words[1], "Set")
                button_remove = PathSettingsDialogButton(words[1], "Remove")
                layout.addWidget(label, idx,0)
                layout.addWidget(file_name_edit, idx+1, 0)
                layout.addWidget(button, idx+1, 1)
                layout.addWidget(button_remove, idx+1, 2)
                button.button_clicked.connect(self.set_path)
                button_remove.button_clicked.connect(self.remove_path)
                idx += 2
        ok_button = QtWidgets.QPushButton("OK")
        ok_button.clicked.connect(self.on_changes_accepted)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.on_cancel)
        ok_cancel_layout = QtWidgets.QHBoxLayout()
        ok_cancel_layout.addWidget(ok_button)
        ok_cancel_layout.addWidget(cancel_button)
        ok_cancel_layout.addStretch(10)
        layout.addLayout(ok_cancel_layout, idx, 0, 1, 2)
        self.setLayout(layout)

    def remove_path(self, name):
        self.paths[name] = ""
        self.path_edits[name].setText("")
        self.settings.remove(name)

    def set_path(self, name):
        if name not in self.files:
            directory = QtWidgets.QFileDialog.getExistingDirectory(caption=name+" directory")
            if directory:
                self.paths[name] = directory
                self.path_edits[name].setText(directory)
        else:
            file_path = QtWidgets.QFileDialog.getOpenFileName(caption=name+" file")
            if file_path:
                self.paths[name] = file_path[0]
                self.path_edits[name].setText(file_path[0])

    def on_changes_accepted(self):
        for program in self.paths.keys():
            current_value = self.settings.get("paths/"+program)
            if current_value != self.paths[program]:
                self.settings.set("paths/"+program, self.paths[program])
        self.close()

    def on_cancel(self):
        self.close()

