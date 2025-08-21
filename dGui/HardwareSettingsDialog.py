from Settings import Settings
from PySide6 import QtWidgets
from PySide6.QtCore import Signal


class HardwareSettingsDialog(QtWidgets.QDialog):
    hardware_settings_updated = Signal()

    def __init__(self):
        super().__init__()
        self.settings = Settings()
        layout = QtWidgets.QGridLayout()
        all_keys = self.settings.keys()
        label = QtWidgets.QLabel("n CPU")
        layout.addWidget(label, 0, 0)
        n_cpu = int(self.settings.get("hardware/nCPU"))
        self.n_cpu_spin_box = QtWidgets.QSpinBox(self)
        self.n_cpu_spin_box.setValue(n_cpu)
        layout.addWidget(self.n_cpu_spin_box, 0, 1)

        memory = self.settings.get("hardware/memory")
        memory_number = int(memory[0:-2])
        memory_units = memory[-2:]

        memory_label = QtWidgets.QLabel("memory")
        layout.addWidget(memory_label, 1, 0)
        self.memory_spin_box = QtWidgets.QSpinBox(self)
        self.memory_spin_box.setValue(memory_number)
        layout.addWidget(self.memory_spin_box, 1, 1)
        self.memory_unit_combo_box = QtWidgets.QComboBox()
        self.memory_unit_combo_box.addItems(["MB", "GB"])
        self.memory_unit_combo_box.setCurrentText(memory_units)
        layout.addWidget(self.memory_unit_combo_box, 1, 2)

        ok_button = QtWidgets.QPushButton("OK")
        ok_button.clicked.connect(self.on_changes_accepted)
        cancel_button = QtWidgets.QPushButton("cancel")
        cancel_button.clicked.connect(self.on_cancel)
        ok_cancel_layout = QtWidgets.QHBoxLayout()
        ok_cancel_layout.addWidget(ok_button)
        ok_cancel_layout.addWidget(cancel_button)
        ok_cancel_layout.addStretch(10)
        layout.addLayout(ok_cancel_layout, 2, 0, 1, 2)
        self.setLayout(layout)

    def on_changes_accepted(self):
        self.settings.set("hardware/nCPU", str(self.n_cpu_spin_box.value()))
        memory = str(self.memory_spin_box.value()) + self.memory_unit_combo_box.currentText()
        self.settings.set("hardware/memory", memory)
        self.hardware_settings_updated.emit()
        self.close()

    def on_cancel(self):
        self.close()
