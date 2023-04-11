
from magicgui import magic_factory
from qtpy.QtWidgets import QWidget

from . import load_path

@magic_factory
def make_widget(file_path: "pathlib.Path" = Path()):
    filename = str(file_path)
    load_path(filename)

