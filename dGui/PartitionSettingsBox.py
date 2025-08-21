from PySide6 import QtWidgets

class PartitionSettingsBox(QtWidgets.QGroupBox):

    def __init__(self, parent = None):
        super().__init__(parent)
        self.layout = QtWidgets.QGridLayout()

        self.exp_har_n_label = QtWidgets.QLabel("exponent")
        self.expHarParameters_doubleSpinBox = QtWidgets.QDoubleSpinBox()
        self.exp_har_widgets = [self.exp_har_n_label, self.expHarParameters_doubleSpinBox]
        self.init_exp_har()

        self.mbis_max_n_steps_label = QtWidgets.QLabel("max number of steps")
        self.mbis_max_n_steps_spin_box =  QtWidgets.QSpinBox()
        self.mbis_convergence_label = QtWidgets.QLabel("convergence threshold")
        self.mbis_convergence_threshold_line_edit = QtWidgets.QLineEdit()
        self.mbis_widgets = [self.mbis_max_n_steps_label, self.mbis_max_n_steps_spin_box, self.mbis_convergence_label,
                             self.mbis_convergence_threshold_line_edit]
        self.init_mbis()

        upper_row_layout = QtWidgets.QHBoxLayout()
        partition_label = QtWidgets.QLabel("partition")
        self.partition_combo_box = QtWidgets.QComboBox()
        self.partition_combo_box.currentTextChanged.connect(self.on_partition_changed)
        self.partition_combo_box.addItems(["Hirshfeld", "exponential Hirshfeld", "MBIS"])
        upper_row_layout.addWidget(partition_label)
        upper_row_layout.addWidget(self.partition_combo_box)
        upper_row_layout.addStretch(10)
        #self.layout.addWidget(partition_label, 0, 0)
        #self.layout.addWidget(self.partition_combo_box, 0, 1, 1, 2)
        self.layout.addLayout(upper_row_layout,0, 0, 1, 2)
        self.layout.setColumnStretch(self.layout.columnCount(), 10)
        self.setLayout(self.layout)

    def init_mbis(self):
        self.mbis_max_n_steps_spin_box.setValue(200)
        self.mbis_max_n_steps_spin_box.setMaximum(100000)
        self.mbis_convergence_threshold_line_edit.setText("0.0001")
        self.layout.addWidget(self.mbis_max_n_steps_label, 1, 0)
        self.layout.addWidget(self.mbis_max_n_steps_spin_box, 1, 1)
        self.layout.addWidget(self.mbis_convergence_label, 2, 0)
        self.layout.addWidget(self.mbis_convergence_threshold_line_edit, 2, 1)

    def init_exp_har(self):
        self.expHarParameters_doubleSpinBox.setValue(1.0)
        self.expHarParameters_doubleSpinBox.setDecimals(3)
        self.expHarParameters_doubleSpinBox.setSingleStep(0.1)
        self.layout.addWidget(self.exp_har_n_label, 1, 0)
        self.layout.addWidget(self.expHarParameters_doubleSpinBox, 1, 1)

    def hide_widgets(self):
        for widget in self.exp_har_widgets:
            widget.hide()
        for widget in self.mbis_widgets:
            widget.hide()

    def on_partition_changed(self, text):
        self.hide_widgets()
        if text == "exponential Hirshfeld":
            for widget in self.exp_har_widgets:
                widget.show()
        if text == "MBIS":
            for widget in self.mbis_widgets:
                widget.show()
