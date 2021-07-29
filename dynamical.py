#currently an archive of an early attempt, not in use.
    def dynamical(self):
        """Dynamical spectra, not implemented"""
        t=tk.Toplevel(self.master,height=600,width=600)
        t.wm_title("Dynamical Spectrum")
        dynamicalframe=tk.Frame(t)
        dynamicalframe.pack(side="top", fill="both",expand=1)
        dynamical=plt.figure()
        ax = dynamical.add_subplot(111)
        ax.tick_params(right= True,top= True,which="both")
        canvas = FigureCanvasTkAgg(dynamical, master=dynamicalframe)
        canvas.get_tk_widget().pack(side="top", fill="both")
        try:
            self.toolbar = NavigationToolbar2TkAgg(canvas, dynamicalframe )
        except:
            self.toolbar = NavigationToolbar2Tk(canvas, dynamicalframe )

        jd=[]
        for i,row in enumerate(self.database):
            self.header=self.database[row]['header']
            try:
                jd.append(self.header['JD'])
            except:
                jd.append(i)
        wavelengthbase=self.database[self.stack[0]]['wavelength']
        x=len(wavelengthbase)

        dynamical.suptitle("PySplot - Date: "+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now(),fontsize=10))
        ax.set_ylabel("JD-%f3"%float(min(jd)))
        # self.xaxislabel()
        tt=np.arange(0, max(jd)-min(jd),.5) #phase steps

        data=np.empty([len(tt),x],dtype=float)
        data.fill(np.nan)
        jd=np.array(jd)-min(jd)
        for i,row in enumerate(self.database):
            flux=self.database[row]['flux']
            wavelength=self.database[row]['wavelength']
            temp_spec=np.empty(x,dtype=float)
            for j,val in enumerate(wavelength):
              temp_spec[j]=np.interp(val,x,flux)
            try:
                data[find_nearest_index(tt,self.jd[i])-1,:]=temp_spec
            except:
                pass
            data[find_nearest_index(tt,self.jd[i]),:]=temp_spec
            try:
                data[find_nearest_index(tt,self.jd[i])+1,:]=temp_spec
            except:
                pass
        # data needs to be resampled so all is on same sampling.
        # interpolate, but only along the time axis
        # for i,junk in enumerate(data[0,:]):
        #     data[:,i]=fill_nan(data[:,i])

        cmap=plt.get_cmap('Spectral')
        # cmap.set_under(color='white')
        # cbaxes = self.fig.add_axes([.88, .32, 0.03, .58])
        # ax.colorbar(cax = cbaxes).ax.tick_params(axis='y', direction='out')  #not sure where defined.

        ax.imshow(data,cmap=cmap, origin='lower',aspect='auto',alpha=1)#, extent=(self.wavelength[1]/u.Angstrom,self.wavelength[-1]/u.Angstrom,tt[0],tt[-1])),interpolation='nearest'
        self.toolbar.update()
        canvas.draw()
