from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, QSize

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.image as img
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.colors as colors

class MainWidget(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(MainWidget, self).__init__(parent)


    def create_mainwidget(self):
        self.create_pysplot()
        self.create_sidepane()
        self.create_output()
        self.parent().mainwidget=QtWidgets.QWidget()
        self.parent().setCentralWidget(self.parent().mainwidget)

        self.grid_layout=QtWidgets.QGridLayout()
        self.grid_layout.addWidget(self.outputbox,0,0)
        self.grid_layout.addWidget(self.pysplotbox,1,0)
        self.grid_layout.addWidget(self.sidepanebox,0,1,2,1)

        self.grid_layout.setColumnStretch(0,4)
        self.grid_layout.setColumnStretch(1,1)
        self.grid_layout.setRowStretch(0,1)
        self.grid_layout.setRowStretch(1,7)


        self.grid_layout.setSpacing(0)
        self.grid_layout.setHorizontalSpacing(0)
        self.grid_layout.setVerticalSpacing(0)

        self.parent().mainwidget.setLayout(self.grid_layout)


    def stackpane(self):
        try:
            self.sidepanebox.close()
        except:
            pass
        self.create_sidepane()
        self.grid_layout.addWidget(self.sidepanebox,0,1,2,1)
        self.parent().mainwidget.setLayout(self.grid_layout)



    def create_sidepane(self):
        self.sidepanebox=QtWidgets.QFrame()
        layout=QtWidgets.QVBoxLayout()
        scroll = QtWidgets.QScrollArea()             # Scroll Area which contains the widgets, set as the centralWidget
        widget = QtWidgets.QWidget()                 # Widget that contains the collection of Vertical Box
        self.vbox = QtWidgets.QVBoxLayout()           # The Vertical Box that contains the Horizontal Boxes of  labels and buttons
        # button1=QtWidgets.QPushButton("Delete Selected")
        # button1.clicked.connect(self.parent().removefromstack)
        self.buildstack()
        # layout.addWidget(button1)
        widget.setLayout(self.vbox)

        #Scroll Area Properties
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        # scroll.setWidgetResizable(False)
        scroll.setWidget(widget)
        layout.addWidget(scroll)
        self.sidepanebox.setLayout(layout)


    def outputupdate(self):
        try:
            self.outputbox.close()
        except:
            pass
        self.create_output()
        self.grid_layout.addWidget(self.outputbox,0,0)
        self.parent().mainwidget.setLayout(self.grid_layout)


    def create_pysplot(self):
        self.pysplotbox = QtWidgets.QFrame()
        self.fig=Figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.tick_params(right= True,top= True,which="both")
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar( self.canvas, self)

        layout = QtWidgets.QVBoxLayout()
        # layout.addWidget(button)
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)

        self.ax.set_title("Single Spectra Mode",fontsize=12)
        # self.logo()
        self.toolbar.update()
        self.canvas.draw()
        # self.pysplotbox=QtWidgets.QWidget()
        self.pysplotbox.setLayout(layout)
        # self.setCentralWidget(widget)
    #
    # def logo(self):
          # '''Display logo on load.  Breaks the display.'''
    #     self.ax.imshow(img.imread('icon.png'))


    def buildstack(self):
        try:
            # specnumber=list(np.arange(len(self.database))+1)
            # specnumber.reverse()
            # i=0
            specforward=list(self.database)
            speclist=specforward.copy()
            speclist.reverse()

            for f in speclist: #the syntax [::-1] reverses the list without modifying it so that  the button list is the same vertical order as the stack plotted spectra.
                b=QtWidgets.QPushButton("%s: %s"%(self.database[f]['stacknumber'],os.path.basename(f)))
                # b=QtQidgets.QPushButton("<font color=%s> %s: %s</font>"%(colors.to_hex(self.database[self.fname]['plotcolor']),t,os.path.basename(s)))
                # b.setText("<font color=%s> %s: %s</font>"%(colors.to_hex(self.database[self.fname]['plotcolor']),t,os.path.basename(s)))

                b.setStyleSheet("Text-align:left")
                # b.setStyleSheet("QPushButton::pressed"
                #                 "{"
                #                 "background-color : red;"
                #                 "}")
                b.clicked.connect(partial(self.stackwindowplot,f))
                # b=QtWidgets.QLabel()
                # b.setBackground(colors.to_hex(self.database[self.fname]['plotcolor']))#QColor(self.database[self.fname]['plotcolor']))
                # b.setColor(colors.to_hex(self.database[self.fname]['plotcolor']))
                # b.setStyleSheet("background-color: white")
                # b.setTextInteractionFlags(Qt.TextSelectableByMouse)
                self.vbox.addWidget(b)
        except:
            pass

    def create_output(self):
        self.outputbox=QtWidgets.QWidget()
        layout=QtWidgets.QVBoxLayout()
        scroll = QtWidgets.QScrollArea()             # Scroll Area which contains the widgets, set as the centralWidget
        widget = QtWidgets.QWidget()                 # Widget that contains the collection of Vertical Box
        vbox = QtWidgets.QVBoxLayout()           # The Vertical Box that contains the Horizontal Boxes of  labels and buttons
        for m in self.parent().message[::-1]: #the minus 1 puts the list inverted because the scrollbox always starts at the top....
            b=QtWidgets.QLabel(m)
            # b.setStyleSheet("background-color: white")
            b.setTextInteractionFlags(Qt.TextSelectableByMouse)
            vbox.addWidget(b)
        widget.setLayout(vbox)

        #Scroll Area Properties
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        scroll.setWidgetResizable(True)
        scroll.setWidget(widget)
        layout.addWidget(scroll)
        self.outputbox.setLayout(layout)
