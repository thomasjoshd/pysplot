from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize
import numpy as np


class ArithWin(QtWidgets.QMainWindow):
    def __init__(self, parent):
        super(ArithWin, self).__init__(parent)
        self.arithsetupUI()
        self.init()

        bar = self.menuBar()
        filemenu = bar.addMenu('File')
        exit=QtWidgets.QAction("Close Arith Window", self)
        exit.setShortcut('Esc')
        exit.triggered.connect(self.quit)
        filemenu.addAction(exit)


    def quit(self):
        self.close()

    def init(self):
        self.slist1=[]
        self.slist2=[]
        self.op=''
        self.sel1=[]
        self.sel2=[]
        self.s=[]
        self.const=0
        self.out=[]

    def reset(self):
        try:
            del self.e1
        except:
            pass
        try:
            del self.e2
        except:
            pass
        try:
            del self.outname
        except:
            pass
        self.arithsetupUI()
        self.init()


    def arithsetupUI(self):
        """Arithmetic"""

        self.setWindowTitle("Spectrum Arithmetic")

        self.top=QtWidgets.QWidget()
        self.setCentralWidget(self.top)

        self.grid=QtWidgets.QGridLayout()
        l1=QtWidgets.QLabel("Operand 1 (Enter Spectrum Number(s) (1,5) or Range (1-5)):")
        h1=QtWidgets.QLabel("Interpolation occurs based on pair of bluest points.")
        self.e1=QtWidgets.QLineEdit()
        l2=QtWidgets.QLabel("Operand 2 (Enter Spectrum Number or Constant):")
        h2=QtWidgets.QLabel("If Operand 2 will be applied to all in the list of operand 1.")
        self.e2=QtWidgets.QLineEdit()
        l3=QtWidgets.QLabel("Operation:")

        l4=QtWidgets.QLabel("Result Name, or suffix:")
        self.outname=QtWidgets.QLineEdit()

        b1=QtWidgets.QPushButton("Go Go Spectra")
        b1.clicked.connect(self.sarith)
        b2=QtWidgets.QPushButton("Reset")
        b2.clicked.connect(self.reset)


        OperationBox = QtWidgets.QComboBox(self)
        # OperationBox.addItem(" ")
        OperationBox.addItem("+ spectrum")
        OperationBox.addItem("- spectrum")
        OperationBox.addItem("/ spectrum")
        OperationBox.addItem("* spectrum")
        OperationBox.addItem("+ constant")
        OperationBox.addItem("- constant")
        OperationBox.addItem("/ constant")
        OperationBox.addItem("* constant")
        OperationBox.addItem("Average Operand 1")
        OperationBox.addItem("Median Operand 1")
        OperationBox.addItem("Add All Operand 1")
        OperationBox.activated[str].connect(self.operation)

        self.grid.addWidget(l1,0,0)
        self.grid.addWidget(self.e1,0,1)
        self.grid.addWidget(h1,1,0)
        self.grid.addWidget(l3,2,0)
        self.grid.addWidget(OperationBox,2,1)
        # self.grid.addWidget(self.e3,1,1)
        self.grid.addWidget(l2,3,0)
        # self.grid.addWidget(O2Box,2,1)
        self.grid.addWidget(self.e2,3,1)
        self.grid.addWidget(h2,4,0)
        self.grid.addWidget(l4,5,0)
        self.grid.addWidget(self.outname,5,1)
        self.grid.addWidget(b1,6,0)
        self.grid.addWidget(b2,6,1)


        # self.grid.addWidget(self.hlist,3,0)
        self.top.setLayout(self.grid)

    def operation(self,text):
        self.op=text
        # print("operand",self.op)

    # def op2(self,text):
    #     print(text,self.e3.text())

    def sarith(self):
        try:
            self.sel1=list(self.e1.text().split(","))
            try:
                self.sel2=list(self.e2.text().split(","))
                self.minilist2()
            except:
                pass
            # self.sel2=self.e2.text()
            self.minilist1()

            # for self.s in self.slist1:
            self.outlist()
            # print(self.out)
            if self.op == '+ spectrum':
                # self.slist2.append(self.parent().stack[int(s)-1])
                self.add()
            elif self.op == '- spectrum':
                # self.slist2.append(self.parent().stack[int(s)-1])
                self.sub()
            elif self.op == '/ spectrum':
                # self.slist2.append(self.parent().stack[int(s)-1])
                self.divspec()
            elif self.op == '* spectrum':
                # self.slist2.append(self.parent().stack[int(s)-1])
                self.multspec()
            elif self.op == '+ constant':
                self.const=float(self.sel2[0])
                self.add()
            elif self.op == '- constant':
                self.const=float(self.sel2[0])
                self.sub()
            elif self.op == '/ constant':
                self.const=float(self.sel2[0])
                self.divide_c()
            elif self.op == '* constant':
                self.const=float(self.sel2[0])
                self.mult_c()
            elif self.op == "Average Operand 1":
                if self.outname.text() != '':
                    self.parent().stackcombine(stack=self.slist1,func='ave',outname=self.outname.text())
                elif self.outname.text() == '':
                    self.parent().stackcombine(stack=self.slist1,func='ave')
            elif self.op == "Median Operand 1":
                if self.outname.text() != '':
                    self.parent().stackcombine(stack=self.slist1,func='median',outname=self.outname.text())
                elif self.outname.text() == '':
                    self.parent().stackcombine(stack=self.slist1,func='median')
            elif self.op == "Add All Operand 1":
                if self.outname.text() != '':
                    self.parent().stackadd(stack=self.slist1,outname=self.outname.text())
                elif self.outname.text() == '':
                    self.parent().stackadd(stack=self.slist1)
            else:
                print("argument error")
            self.close()
        except:
            print("Problem parsing")
            print(self.op,self.sel1,self.sel2,self.outname.text(),self.out,self.slist1,self.slist2)
            self.reset()

    def outlist(self):
        self.out=[]
        if self.outname.text() == '':
            suffix='1'
        else:
            suffix=self.outname.text()
        for val in self.slist1:
            self.out.append(val+suffix)

    def divide_c(self):
        for i,val in enumerate(self.slist1):
            print(self.out[i])
            self.parent().divide_const(spec=val,divisor=self.const,outname=self.out[i])

    def mult_c(self):
        for i,val in enumerate(self.slist1):
            self.parent().mult_const(spec=val,multiplier=self.const,outname=self.out[i])

    def divspec(self):
        #divide spectra by spectrum
        for i,val in enumerate(self.slist1):
                a=[]
                a.append(val)
                a.append(self.slist2[0])
                self.parent().dividespec(stack=a,outname=self.out[i])

    def multspec(self):
        #divide spectra by spectrum
        for i,val in enumerate(self.slist1):
                a=[]
                a.append(val)
                a.append(self.slist2[0])
                self.parent().multiplyspec(stack=a,outname=self.out[i])


    def add(self):
        #Add all spectra in list together, or add a constant
        if abs(self.const) > 0: #add constant to spectra in list
            for i,a in enumerate(self.slist1):
                self.parent().stackadd_const(stack=[a],const=self.const,outname=self.out[i])
        else: #add all spectra together
            # if len(self.slist2) > 1:
            #     a=[]
            #     for val in self.slist1:
            #         a.append(val)
            #     for val in self.slist2:
            #         a.append(val)
            #     if self.outname.text() != '' and float(self.const) == 0:
            #         self.parent().stackadd(stack=a,outname=self.outname.text())
            #     elif self.outname.text() == '' and float(self.const) == 0:
            #         self.parent().stackadd(stack=a)
            # else: # add operand 2 spectrum to all the spectra in the first list.
            for i,val in enumerate(self.slist1):
                a=[]
                a.append(val)
                a.append(self.slist2[0])
                self.parent().stackadd(stack=a,outname=self.out[i])


    def sub(self):
        if abs(self.const) > 0: #subtract a constant from each spectrum in list
            for i,a in enumerate(self.slist1):
                self.parent().stackadd_const(stack=[a],const=-self.const,outname=self.out[i])
        else: #subtract one spectrum from another
            for i,val in enumerate(self.slist1):
                a=[]
                a.append(val)
                a.append(self.slist2[0])
                self.parent().stacksubtract(stack=a,outname=self.out[i])


    def minilist1(self):
        # print("minilist1")
        for s in self.sel1:
            if '-' in s:
                x,y=s.split("-")
                x=int(x)
                y=int(y)
                # print(type(x))
                srange=np.arange(x,y+1,1,dtype=int)
                # srange=np.arange(1,6)
                # print('hello')
                # print(srange)
                for k in srange:
                    # print(k)
                    self.slist1.append(self.parent().stack[int(k)-1])
            elif s == '0':
                pass
            else:
                self.slist1.append(self.parent().stack[int(s)-1])
    #
    def minilist2(self):
        # print("minilist2")
        for s in self.sel2:
            if '-' in s:
                x,y=s.split("-")
                x=int(x)
                y=int(y)
                # print(type(x))
                srange=np.arange(x,y+1,1,dtype=int)
                # srange=np.arange(1,6)
                # print(srange)
                for k in srange:
                    # print(k)
                    self.slist2.append(self.parent().stack[int(k)-1])
            elif s == '0':
                pass
            else:
                self.slist2.append(self.parent().stack[int(s)-1])
