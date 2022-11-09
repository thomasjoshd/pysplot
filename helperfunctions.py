import numpy as np

def find_nearest_index(A,value):
    """Finds the closest data points to the input array."""
    return (np.abs(np.asarray(A)-value)).argmin()

def julian(date):
    """Returns the Julian Date of a Gregorian Calendar Date."""
    year,month,day,UTh,UTm,UTs=date
    a = (14 - month)//12
    y = year + 4800 - a
    m = month + 12*a - 3
    UT=(UTh-12)/24+UTm/1440+UTs/86400
    return day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045 +UT

    # def df1(self,x):
    #     return np.interp(x,self.xm1,self.diff1)
    # def df2(self,x):
    #     return np.interp(x,self.xm2,self.diff2)
    # def hf(self):
    #     return self.height

def snr(flux):
    try:
        mu=np.average(flux)
        sigma=np.std(flux)
        return mu/sigma
    except:
        print("Excpetion in helperfunction.snr")
