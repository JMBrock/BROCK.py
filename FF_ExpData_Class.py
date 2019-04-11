class Experimental_Dynamic_Data:
    """
        This class is created to read in, fit, and derive various dynamic data from an external
        experimental data source
    """

    def __init__(self,datafile):
        self.dfile = datafile

    def startup(self):
        import numpy as np
        # Load in data
        print("-- Reading in '%s' File" % self.dfile)
        self.ExpData = np.loadtxt(self.dfile)

    def startup_MultiBody(self):
        import numpy as np
        import ProgressBar

        # Load in data
        print("\n-- Reading in '%s' File" % self.dfile)

        with open(self.dfile, 'r', errors='replace') as f:
            lines_ = f.readlines()
        lines = lines_[3:]

        self.Ndata = len(lines)
        self.ExpData = np.zeros((self.Ndata,7),dtype=np.float64)
        for i in range(self.Ndata):
            self.ExpData[i][:] = lines[i].split()
            precent = float(i+1)/float(self.Ndata)*100.0
            ProgressBar.ProgressBar(i+1,self.Ndata,precent)


    def functional_fit(self,XYZFitType,PitchFitType,YawFitType,NFit):
        from scipy.interpolate import interp1d
        import numpy as np
        """
            Fit the incoming data
        """
        self.NFit = NFit
        self.XYZFitType = XYZFitType
        self.PitchFitType = PitchFitType
        self.YawFitType = YawFitType


        self.Fit_X = interp1d( self.Time, self.X_loc, kind=XYZFitType )
        self.Fit_Y = interp1d( self.Time, self.Y_loc, kind=XYZFitType )
        self.Fit_Z = interp1d( self.Time, self.Z_loc, kind=XYZFitType )
        self.Fit_Yaw = interp1d( self.Time, self.Yaw, kind=YawFitType )
        self.Fit_Pitch = interp1d( self.Time, self.Pitch,kind=PitchFitType )

        self.Time_Fit = np.linspace(self.Time[0], self.Time[self.Ndata-1], num=self.NFit, endpoint=True)

        self.X_loc_Fit = self.Fit_X(self.Time_Fit)
        self.Y_loc_Fit = self.Fit_Y(self.Time_Fit)
        self.Z_loc_Fit = self.Fit_Z(self.Time_Fit)
        self.Yaw_Fit = self.Fit_Yaw(self.Time_Fit)
        self.Pitch_Fit = self.Fit_Pitch(self.Time_Fit)

    def estimate_rates(self):
        import numpy as np
        """
            Estimate tip-off rates
        """

        rad2deg = 180.0/np.pi

        deltaT = self.Time_Fit[1] - self.Time_Fit[0]
        delta_X = self.X_loc_Fit[1] - self.X_loc_Fit[0]
        delta_Y = self.Y_loc_Fit[1] - self.Y_loc_Fit[0]
        delta_Z = self.Z_loc_Fit[1] - self.Z_loc_Fit[0]
        delta_Yaw = self.Yaw_Fit[1] - self.Yaw_Fit[0]
        delta_Pitch = self.Pitch_Fit[1] - self.Pitch_Fit[0]

        print( "-- Estimated Fitted Estimate rate for pitch: %f" % ( (delta_Pitch/deltaT)*rad2deg) )
        print( "-- Estimated Fitted Estimate rate for yaw:   %f" % ( (delta_Yaw/deltaT)*rad2deg) )
        print( "-- Linear x-velocity: %f" % ( (delta_X/deltaT) ) )
        print( "-- Linear y-velocity: %f" % ( (delta_Y/deltaT) ) )
        print( "-- Linear z-velocity: %f" % ( (delta_Z/deltaT) ) )

    def Extract_Raw_Data(self,Time_id,X_id,Y_id,Z_id,Pitch_id,Yaw_id,Roll_id):
        """
            Read in experimental data. This routine requires an input for the time,
            x,y,z,pitch,yaw, and roll index IDs.
        """

        # zero out values at the start of the run
        self.TimeRaw = self.ExpData[:,Time_id]
        self.X_locRaw = self.ExpData[:,X_id]
        self.Y_locRaw = self.ExpData[:,Y_id]
        self.Z_locRaw = self.ExpData[:,Z_id]
        self.Time  = self.ExpData[:,Time_id] - self.ExpData[0,Time_id]
        self.X_loc = self.ExpData[:,X_id] - self.ExpData[0,X_id]
        self.Y_loc = self.ExpData[:,Y_id] - self.ExpData[0,Y_id]
        self.Z_loc = self.ExpData[:,Z_id] - self.ExpData[0,Z_id]
        self.Pitch = self.ExpData[:,Pitch_id]
        self.Yaw   = self.ExpData[:,Yaw_id]
        self.Roll  = self.ExpData[:,Roll_id]
        self.Ndata = len(self.Time)

    def Wash_Data(self):
        import numpy as np
        """
           Clean up indexing for cases with bad data (we will skip them)
           This needs some work since there are gonna be some stations were
           one camera has good data and the other doesn't...
        """
        
        self.pindex=[]
        self.yindex=[]
        for i in range(self.Ndata):
            if self.Pitch[i]==-0.0:
                self.pindex.append(i)
            if self.Yaw[i]==-0.0:
                self.yindex.append(i)

        self.Pitch = np.delete(pitch,pindex)
        self.Time = np.delete(t,pindex)
        self.Yaw = np.delete(yaw,yindex)

    def Compute_X_VoT(self):
        """
            Compute X-V0T for frame displacement values in the trajectory file
        """

        # Estimate V0T if not given
        deltaT = self.Time_Fit[1] - self.Time_Fit[0]
        delta_X = self.X_loc_Fit[1] - self.X_loc_Fit[0]

        V0 = delta_X/deltaT

        self.X_VoT = 0.0*self.X_loc_Fit
        VoT = 0.0
        for i in range(self.NFit):
            self.X_VoT[i] =  (self.X_loc_Fit[i] - V0*self.Time_Fit[i])

        # self.X_VoT[0] = 0.0

    def Write_TrajectoryFile(self,fname,header,Type):


        # Create a single body trajectory file but with frame discplacement
        print("Creating Trajectory File %s" % (fname))

        f1 = open(str(fname),"wt")
        if header:
            f1.write('#Variables="time",')
            f1.write('"xloc",')
            f1.write('"yloc",')
            f1.write('"zloc",')
            f1.write('"theta",')
            f1.write('"psi",')
            f1.write('"phi"\n')

        for i in range(self.NFit):
            f1.write("%.9f      " %self.Time_Fit[i])
            if Type == 'Lagrange':
                f1.write("%.9f      " %self.X_VoT[i])
            if Type == 'Euler':
                f1.write("%.9f      " %self.X_loc_Fit[i])
            f1.write("%.9f      " %self.Y_loc_Fit[i])
            f1.write("%.9f      " %self.Z_loc_Fit[i])
            f1.write("%.9f      " %self.Pitch_Fit[i])
            f1.write("%.9f      " %self.Yaw_Fit[i])
            f1.write("%.9f    \n" %0.0) # Keeping this zero for now

        f1.close()

    def TikzData(self,fname):
        import numpy as np

        print("-- Creating Tikz file '%s'" % fname)
        f = open(fname,"wt")

        sim_alpha = self.Pitch
        sim_beta = self.Yaw
        self.TotalAlpha = np.arccos(np.cos(self.Pitch*np.pi/180.0)*np.cos(self.Yaw*np.pi/180.0))*180.0/np.pi


        delta_X = self.X_loc[1] - self.X_loc[0]
        deltaT = self.Time[1] - self.Time[0]
        Vmag = delta_X/deltaT
        Vmag = 424.5

        self.VoT = self.X_loc - Vmag*self.Time

        for j in range(self.Ndata):
            f.write("%.9f      " %(self.Time[j])) # time(s)
            f.write("%.9f      " %(self.VoT[j]))  #  x
            f.write("%.9f      " %(self.Y_loc[j]))  #  y
            f.write("%.9f      " %(self.Z_loc[j]))  #  z
            f.write("%.9f      " %(self.Pitch[j])) #  pitch
            f.write("%.9f      " %(self.Yaw[j]))  #  yaw
            f.write("%.9f      " %(0.0))  #  roll
            f.write("%.9f      " %(1.0))  #  error
            f.write("%.9f\n" %(self.TotalAlpha[j]))    #  total alpha
        f.close()