from enum import IntEnum


class AtomDisplayStyle(IntEnum):
    line = 0
    stick = 1
    ball_and_line = 2
    ball_and_stick = 3
    cpk = 4
    ellipsoid_and_line = 5
    ellipsoid_and_stick = 6

    def show_text(self):
        return self.name.replace("_", " ")

    @staticmethod
    def default():
        return AtomDisplayStyle.ellipsoid_and_stick

    @staticmethod
    def from_show_text(text):
        return AtomDisplayStyle[text.replace(" ", "_")]

    @staticmethod
    def get_show_text_list():
        return [key.replace("_", " ") for key in AtomDisplayStyle.__members__.keys()]


