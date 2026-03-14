import os
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize

class HeaderWin(QtWidgets.QMainWindow):
    def __init__(self, parent):
        super(HeaderWin, self).__init__(parent)
        self.setupUI()

        bar = self.menuBar()
        filemenu = bar.addMenu('File')
        exit=QtWidgets.QAction("Close Header Window", self)
        exit.setShortcut('Esc')
        exit.triggered.connect(self.quit)
        filemenu.addAction(exit)


    def quit(self):
        self.close()


    def setupUI(self):
        """Display the FITS header as loaded from the file."""


        self.setWindowTitle("FITS Header")
        # self.setStyleSheet("background-color: white")
        self.resize(800,800)

        self.top=QtWidgets.QWidget()
        self.setCentralWidget(self.top)

        self.grid=QtWidgets.QGridLayout()
        l0=QtWidgets.QLabel("Spectrum:")
        l01=QtWidgets.QLabel(os.path.basename(self.parent().fname))
        l1=QtWidgets.QLabel("Keyword:")
        self.e1=QtWidgets.QLineEdit()
        l2=QtWidgets.QLabel("Value:")
        self.e2=QtWidgets.QLineEdit()
        l3=QtWidgets.QLabel("Comment:")
        self.e3=QtWidgets.QLineEdit()

        b1=QtWidgets.QPushButton("Set Value")
        b1.clicked.connect(self.header_set)
        b2=QtWidgets.QPushButton("Set New Keyword")
        b2.clicked.connect(self.newheader_set)
        b3=QtWidgets.QPushButton("Delete Keyword")
        b3.clicked.connect(self.delkey)

        self.hbox()

        self.grid.addWidget(l0,0,0,1,1)
        self.grid.addWidget(l01,0,1,1,1)
        self.grid.addWidget(l1,1,0,1,1)
        self.grid.addWidget(self.e1,1,1,1,1)
        self.grid.addWidget(l2,2,0,1,1)
        self.grid.addWidget(self.e2,2,1,1,1)
        self.grid.addWidget(l2,2,0,1,1)
        self.grid.addWidget(self.e2,2,1,1,1)
        self.grid.addWidget(l3,3,0,1,1)
        self.grid.addWidget(self.e3,3,1,1,1)
        self.grid.addWidget(b1,4,0,1,1)
        self.grid.addWidget(b2,4,1,1,1)
        self.grid.addWidget(b3,4,2,1,1)
        self.grid.addWidget(self.hlist,5,0,1,3)
        self.top.setLayout(self.grid)


        # self.grid.setColumnStretch(0,1)
        # self.grid_layout.setColumnStretch(1,1)
        # self.grid.setRowStretch(4,1)
        # self.grid_layout.setRowStretch(1,7)


    def hbox(self):
        # self.headerbox=QtWidgets.QWidget()
        self.hlist=QtWidgets.QListWidget()
        # print(self.parent().header.keys)
        for key in self.parent().header.keys():
            self.hlist.addItem("%s = %s / %s"%(key,self.parent().header[key],self.parent().header.comments[key]))
        # for keys in self.parent().header.tostring(sep=',').split(','):
        #     self.hlist.addItem("%s "%(keys))

        # print(self.parent().header.tostring(sep=',').split(','))

        self.hlist.itemSelectionChanged.connect(self.grab_keyword)

        self.grid.addWidget(self.hlist,5,0,1,3)
        self.top.setLayout(self.grid)

    def header_set(self):
        """Write the keyword/value to the header stored in the dictionary."""
        try:
            self.parent().header[self.key]=(self.e2.text(),self.e3.text())
            self.parent().Spectra.updatespectrum()
        except:
            print('Exception in HeaderWin.header_set')
        try:
            self.hlist.close()
        except:
            print('Exception in HeaderWin.header_set() while trying to close hlist')
        self.hbox()

    def newheader_set(self):
        """set a new header value"""
        self.parent().header[self.e1.text()]=(self.e2.text(),self.e3.text())
        self.parent().Spectra.updatespectrum()
        try:
            self.hlist.close()
        except:
            pass
        self.hbox()



    def grab_keyword(self):
        """grabs a header keyword, and set it's values to the display boxes"""
        self.key=self.hlist.selectedItems()[0].text().split('=')[0].strip()
        try:
            self.e1.setText(self.key)
        except:
            print('Exception in HeaderWin.grab_keyword grabbing keyword')
        try:
            self.e2.setText(str(self.parent().header[self.key]))#self.e2.setText(line[1].split('/')[0].strip())
        except:
            pass
        try:
            self.e3.setText(str(self.parent().header.comments[self.key]))#self.e2.setText(line[1].split('/')[0].strip())
        except:
            pass


    def delkey(self):
        try:
            self.parent().header.remove(self.key)
        except:
            print('failed to delete keyword from header')
        try:
            self.hlist.close()
        except:
            pass
        self.hbox()
