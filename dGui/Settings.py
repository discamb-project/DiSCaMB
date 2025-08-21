from PySide6.QtCore import QSettings, QCoreApplication


def make_default_settings():
    settings = {
        "label/color": [0, 255, 255],
        "label/font": 'Arial',
        "label/size": 16,
        "paths/orca": None,
        "paths/gaussian": None,
        "paths/MATTS": None,
        "hardware/nCPU": 1,
        "hardware/memory": "2GB"
    }
    return settings


class Settings:
    def __init__(self):
        QCoreApplication.setOrganizationName("DiSCaMB")
        QCoreApplication.setOrganizationDomain("discamb.org")
        QCoreApplication.setApplicationName("dGui")
        self.settings = QSettings()
        self.__combine_settings()

    def __combine_settings(self):
        settings_keys = self.settings.allKeys()
        default_settings = make_default_settings()
        default_settings_keys = default_settings.keys()
        for key in default_settings_keys:
            if key not in settings_keys:
                self.settings.setValue(key, default_settings[key])

    def keys(self):
        return self.settings.allKeys()

    def set(self, key, value):
        self.settings.setValue(key, value)

    def get(self, key):
        return self.settings.value(key)

    def remove(self, key):
        self.settings.remove(key)

    def has_key(self, key):
        return key in self.keys()

