class Simulation_Dynamic_Data:
    """
    This class reads in simulation dynamic data and does work on the data set 
    such as; Compute AeroCoefficients, reduce/coarsen dataset for tikz plotting,
    fitting functions to back out frequency information, etc. 
    """

    def __init__(self,dfile):
        self.dfile = dfile

    def readData(self):
        import numpy as np

        # Load in data
        print("-- Reading in '%s' File" % self.dfile)
        self.Data = np.loadtxt(self.dfile)


    def readBadData(self):
        from RepairRBFFile import brokenRBFile
        # Load in data
        print("-- Reading in '%s' File" % self.dfile)
        self.Data = brokenRBFile(self.dfile)
        self.Data.file_len()
        print("-- Repairing file and writing new '%s' File" % (self.dfile+'_Fixed'))
        self.Data.repair_and_write(self.dfile+'_Fixed',8)
        self.dfile = self.dfile+'_Fixed'

    def ExtractData_SingleBody(self,Time_id,X_id,Y_id,Z_id,Pitch_id,Yaw_id,Roll_id):
        """
            Read in data. This routine requires an input for the time,
            x,y,z,pitch,yaw, and roll index IDs.
        """

        # zero out values at the start of the run
        self.Time  = self.Data[:,Time_id] - self.Data[0,Time_id]
        self.X_loc = self.Data[:,X_id] - self.Data[0,X_id]
        self.Y_loc = self.Data[:,Y_id] - self.Data[0,Y_id]
        self.Z_loc = self.Data[:,Z_id] - self.Data[0,Z_id]
        self.Pitch = self.Data[:,Pitch_id]
        self.Yaw   = self.Data[:,Yaw_id]
        self.Roll  = self.Data[:,Roll_id]
        self.Ndata = len(self.Time)

    def CoarsenData(self,factor):
        import numpy as np
        """
            Downselect data to make it easier for plotting
        """
        self.TimeC = self.Time[::factor]
        self.X_locC = self.X_loc[::factor]
        self.Y_locC = self.Y_loc[::factor]
        self.Z_locC = self.Z_loc[::factor]
        self.PitchC = self.Pitch[::factor]
        self.YawC   = self.Yaw[::factor]
        self.RollC  = self.Roll[::factor]
        self.NdataC = len(self.TimeC)

    def Compute_X_VoT(self):
        """
            Compute X-V0T for frame displacement values in the trajectory file
        """

        # Estimate V0T if not given
        deltaT = self.Time[1] - self.Time[0]
        delta_X = self.X_loc[1] - self.X_loc[0]

        V0 = delta_X/deltaT

        self.X_VoT = 0.0*self.X_loc
        VoT = 0.0
        for i in range(len(self.Time)):
            self.X_VoT[i] =  (self.X_loc[i] - V0*self.Time[i])

    def TikzData(self,fname,factor):
        import numpy as np

        self.CoarseData = self.CoarsenData(factor)
        print("-- Creating Tikz file '%s' \n" % fname)
        f = open(fname,"wt")

        sim_alpha = self.PitchC
        sim_beta = self.YawC
        self.TotalAlpha = np.arccos(np.cos(self.PitchC*np.pi/180.0)*np.cos(self.YawC*np.pi/180.0))*180.0/np.pi

        delta_X = self.X_locC[1] - self.X_locC[0]
        deltaT = self.TimeC[1] - self.TimeC[0]
        Vmag = delta_X/deltaT

        # self.VoT = self.X_locC - Vmag*self.TimeC
        self.VoT = self.X_locC #- Vmag*self.TimeC

        for j in range(self.NdataC):
            f.write("%.9f      " %(self.TimeC[j])) # time(s)
            f.write("%.9f      " %(self.PitchC[j])) # simulation pitch
            f.write("%.9f      " %(self.YawC[j]))  # simulation yaw
            f.write("%.9f      " %(-self.VoT[j]))  # simulation x-VoT 
            f.write("%.9f      " %(self.TotalAlpha[j]))    # simulation total alpha
            f.write("%.9f \n   " %(self.RollC[j]))  # simulation roll
        f.close()