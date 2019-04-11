import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.signal import blackman
from scipy import integrate

class myfft:
    """
    This class is created for the computation of Fast Fourier Transforms 
    on signals outputted by US3D (Or any signal in the future).
    """
    def __init__(self,dfile,vfile):
        self.dfile = dfile
        self.vfile = vfile

    def file_len(self):
        """
        This reads the length of a file
        """
        with open(self.dfile) as f:
            for i, l in enumerate(f):
                pass
        self.nlines = i + 1


    def file_vars(self):
        self.variables = []
        f = open(self.vfile,'r')
        line = f.readline().strip().strip('variables=')
        self.variables = line.split(',')
        f.close

    def file_read(self):
        """
        This reads the data of a file
        """
        self.data = np.zeros((self.nlines-1,len(self.variables)))
        with open(self.dfile,'r') as f:
            lines_ = f.readlines()
        lines = lines_[1:]
        for i in range(len(lines)):
            vars = lines[i].split()
            for j in range(len(vars)):
                v = vars[j].strip()
                try:
                    self.data[i][j] = float(v)
                except ValueError:
                    pass
                    self.data[i][j] = 0.0

if __name__ == "__main__":

#   pubdir = "/home/jmbrock2/PUBLIC"
#   dynamicdir = "SIAD-R_DYNAMIC_RUN/Shot_2643"
    pubdir = "./"
    dynamicdir = ""
    dfile = "dynamics_log.dat"
    vfile = "dynamics_vars.txt"

    #"v-inlet01"   0.244568525525 294.35 294.35 1131.55843604 !Shot 2643
    rho = 0.244568525525
    vmag = 1043.0

    # Gas properties
    universal_gas_constant = 8.3145
    molar_mass_air = 0.0289645
    air_gas_constant = universal_gas_constant/molar_mass_air
    gamma = 1.4

    # Normalization constants
    Diam = 1.4003*2.54
    Rad = Diam/2.0
    Area = np.pi*Rad**2

    # Freestream dynamic pressure q=1/2*rho*V^2
    dynp = 0.5*rho*vmag**2

    dname = os.path.join(pubdir,dynamicdir,dfile)
    vname = os.path.join(pubdir,dynamicdir,vfile)
    print dname

    dyndata = myfft(dname,vname)
    dyndata.file_len()
    dyndata.file_vars()
    dyndata.file_read()
    print " "
    print 'Total number of lines', dyndata.nlines - 1
    print 'Number of variables: ',len(dyndata.variables)
    print 'Variables = ', dyndata.variables[:]
    # Should be variables=time,Fx,Fy,Fz,Mx,My,Mz,Vx,Vy,Vz,a1,a2,a3,w1,w2,w3

    time = np.zeros((dyndata.nlines-1),dtype=np.float64)

    Fx = np.zeros((dyndata.nlines-1),dtype=np.float64)
    Fy = np.zeros((dyndata.nlines-1),dtype=np.float64)
    Fz = np.zeros((dyndata.nlines-1),dtype=np.float64)

    Vx = np.zeros((dyndata.nlines-1),dtype=np.float64)
    Vy = np.zeros((dyndata.nlines-1),dtype=np.float64)
    Vz = np.zeros((dyndata.nlines-1),dtype=np.float64)

    Mx = np.zeros((dyndata.nlines-1),dtype=np.float64)
    My = np.zeros((dyndata.nlines-1),dtype=np.float64)
    Mz = np.zeros((dyndata.nlines-1),dtype=np.float64)

  
    a1 = np.zeros((dyndata.nlines-1),dtype=np.float64)
    a2 = np.zeros((dyndata.nlines-1),dtype=np.float64)
    a3 = np.zeros((dyndata.nlines-1),dtype=np.float64)

    w1 = np.zeros((dyndata.nlines-1),dtype=np.float64)
    w2 = np.zeros((dyndata.nlines-1),dtype=np.float64)
    w3 = np.zeros((dyndata.nlines-1),dtype=np.float64)

    for j in range(dyndata.nlines-1):

        time[j] = float(dyndata.data[j][dyndata.variables.index("time")])

        Fx[j] = float(dyndata.data[j][dyndata.variables.index("Fx")])
        Fy[j] = float(dyndata.data[j][dyndata.variables.index("Fy")])
        Fz[j] = float(dyndata.data[j][dyndata.variables.index("Fz")])

        Mx[j] = float(dyndata.data[j][dyndata.variables.index("Mx")])
        My[j] = float(dyndata.data[j][dyndata.variables.index("My")])
        Mz[j] = float(dyndata.data[j][dyndata.variables.index("Mz")])

        Vx[j] = float(dyndata.data[j][dyndata.variables.index("Vx")])
        Vy[j] = float(dyndata.data[j][dyndata.variables.index("Vy")])
        Vz[j] = float(dyndata.data[j][dyndata.variables.index("Vz")])

        a1[j] = float(dyndata.data[j][dyndata.variables.index("a1")])
        a2[j] = float(dyndata.data[j][dyndata.variables.index("a2")])
        a3[j] = float(dyndata.data[j][dyndata.variables.index("a3")])

        w1[j] = float(dyndata.data[j][dyndata.variables.index("w1")])
        w2[j] = float(dyndata.data[j][dyndata.variables.index("w2")])
        w3[j] = float(dyndata.data[j][dyndata.variables.index("w3")])
    
    N_data = int(dyndata.nlines-1)
    delta_t = float(dyndata.data[N_data-1][0]-dyndata.data[N_data-2][0])
    Samp_Freq = float(1.0/(delta_t))
    Start_Time = float(dyndata.data[0,0])
    Final_Time = float(dyndata.data[N_data-1,0])
    Run_Time = Final_Time - Start_Time
    Time = np.linspace(0.0,Run_Time,N_data)
    print " "
    print 'Initial velocity: %.5e (m/s)' %vmag
    print 'Total Time of Statistics: %.5e (s)' % Run_Time
    print 'Delta T: %.5e (s)' % delta_t
    print 'Sample Frequency: %3.2f (kHz)' % float(Samp_Freq/1000.0)

    xdist = np.zeros((dyndata.nlines-1),dtype=np.float64)
    ydist = np.zeros((dyndata.nlines-1),dtype=np.float64)
    zdist = np.zeros((dyndata.nlines-1),dtype=np.float64)
    VoT   = np.zeros((dyndata.nlines-1),dtype=np.float64)
    # Compute Cross range distance
    xdisto = 0.0
    ydisto = 0.0
    zdisto = 0.0
    timeo = time[0]
    VoTo = 0.0
    for j in range(N_data):
        xdist[j] = xdisto + (vmag-Vx[j])*(time[j]-timeo)
        ydist[j] = ydisto - Vy[j]*(time[j]-timeo)
        zdist[j] = zdisto - Vz[j]*(time[j]-timeo)
        VoT[j] =  VoTo + vmag*(time[j]-timeo)
        xdisto = xdist[j]
        ydisto = ydist[j]
        zdisto = zdist[j]
        timeo = time[j]
        VoTo = VoT[j]

    sim_alpha = a1
    sim_beta = a2
    sim_roll = a3
    sim_xvot = xdist-VoT
    sim_at = np.arccos(np.cos(a1*np.pi/180.0)*np.cos(a2*np.pi/180.0))*180.0/np.pi

    f = open("2623_dyndata.dat","wt")
    k = 0
    for j in range(N_data):
        xdist = xdisto + (vmag-Vx[j])*(time[j]-timeo)
        ydist = ydisto - Vy[j]*(time[j]-timeo)
        zdist = zdisto - Vz[j]*(time[j]-timeo)
        VoT =  VoTo + vmag*(time[j]-timeo)
        xdisto = xdist
        ydisto = ydist
        zdisto = zdist
        timeo = time[j]
        VoTo = VoT
        runtime = time[j]-time[0]
        if(j==k*(N_data/5000)):
           k = k+1
           f.write("%.9f      " %(runtime)) # time(s)
           f.write("%.9f      " %(sim_alpha[j])) # simulation pitch
           f.write("%.9f      " %(sim_beta[j]))  # simulation yaw
           f.write("%.9f      " %(sim_xvot[j]))  # simulation x-VoT 
           f.write("%.9f      " %(sim_at[j]))    # simulation total alpha
           f.write("%.9f \n   " %(sim_roll[j]))  # simulation yaw
