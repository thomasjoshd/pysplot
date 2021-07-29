# import os
import sys

from PyQt5 import QtWidgets

from MainWin import MainWin

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    # sys.excepthook = traceback.print_exception()
    # creating main window
    mw = MainWin()
    mw.show()
    sys.exit(app.exec())
