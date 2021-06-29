import os
from PyQt5 import QtWidgets
import datetime

class LogFiles(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(LogFiles, self).__init__(parent)

    def captainslog(self):
        """Save the output from measuremetns to a csv"""
        try:
            path=os.path.dirname(self.parent().fname)
            time='{0:%Y%m%d.%H%M%S}'.format(datetime.datetime.now())
            basename=os.path.basename("pysplot-%s.csv"%time)
            savename=os.path.join(path,basename)
            self.logname=savename
            self.parent().message.append("Log being written to: %s"%self.logname)
            self.parent().outputupdate()
            self.starlog=open(savename,'w')
            self.starlog.write('PySplot Log %s'%time)
            self.starlog.write('\n')
        except:
            print('Exception occured in LogFiles.captainslog')

    def endoflog(self):
        """close the log file"""
        try:
            self.starlog.close()
            # print(self.logname)
            self.parent().message.append("Log written to: %s"%self.logname)
            self.parent().outputupdate()
            del self.starlog
        except:
            print('Exception in LogFiles.endoflog')

    def checklog(self):
        """checks to see if a log is open."""
        try:
            self.starlog
        except:
            self.captainslog()
            # print('Except, trying creating log, opening Captain's Log')


    def write(self,text):
        try:
            # print(text)
            self.starlog.write(text)
            self.starlog.write('\n')
        except:
            print('Error in LogFiles.write')
