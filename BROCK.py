"""
Python code to read experimental data and create Trajectory file and initial rate estimates. 
   The rates can be used in the older single body version of the FFCFD code while the trajectory.dat 
   file can be used in the newer multibody RBF FFCFD code. Additionally, this suit can read the output
   from either versions of the code and create plots or reduced data for tikz plotting for presentations 
   and papers.

Author: JMB
Date Last Edit: 04-01-2019 minor cleanup
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from FF_ExpData_Class import Experimental_Dynamic_Data
from FF_SimData_Class import Simulation_Dynamic_Data


def main(argv):

    print("%=============================================================%")
    print("|                                                             |")
    print("| Welcome to Ballistic Range on Computer Kernel (BROCK) tool! |")
    print("|                Developed by: Joseph M. Brock                |")
    print("|                                                             |")
    print("%=============================================================%")

    try:
        import argparse
    except ImportError:
        sys.exit("Unable to import argparse!")

    #Read in command line arguments

    parser = argparse.ArgumentParser()

    parser.add_argument('-ExpData', metavar='', type=str, nargs='+',
                        help='input expreimental data file\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-dfile="shot_data.dat"')
    parser.add_argument('-FitData', metavar='', type=str,
                        help='Option to tell program we want to fit data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-FitData')
    parser.add_argument('-CreateTrajectory', action='store_true',
                        help='Option to export trajectory file(s).\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-CreateTrajectory'
                        '\n NOTE: This will create two files; trajectory.dat and trajectory.dat2, '
                        '\n       which will have the true x and x-VoT trajectories, respectively ')
    parser.add_argument('-TrajFile', metavar='', type=str,
                        help='output trajectory data file to be used in US3D\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-TrajFile="trajectory.dat"')
    parser.add_argument('-PitchFitType', metavar='', type=str,
                        help='type of fit function used to fit data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-PitchFitType="cubic/quadratic/slinear"')
    parser.add_argument('-YawFitType', metavar='', type=str,
                        help='type of fit function used to fit data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-YawFitType="cubic/quadratic/slinear"')
    parser.add_argument('-XYZFitType', metavar='', type=str,
                        help='type of fit function used to fit data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-XYZFitType="cubic/quadratic/slinear"')
    parser.add_argument('-NFit', type=int,
                        help='Number of points to use for function fit\n')
    parser.add_argument('-PlotExp', action='store_true',
                        help='Show experimental data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-PlotExp')
    parser.add_argument('-PlotFits', action='store_true',
                        help='Show fitted data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-PlotFits')
    parser.add_argument('-mm2m',action='store_true')
    parser.add_argument('-rad2deg',action='store_true')
    parser.add_argument('-SimData', metavar='', type=str,
                        help='input dynamic data file\n'
                        'example usage: python Get_Rates_Make_Trajectory '
                        ' -SimData=dynamic_log.dat')
    parser.add_argument('-RepairSimData',action='store_true',
                        help='Dynamic datafile from simulation has weird output on multiple lines.'
                        ' We can fix this...')
    parser.add_argument('-PlotSim', action='store_true',
                        help='Show simulation data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-PlotSim')
    parser.add_argument('-PlotCompare', action='store_true',
                        help='Show simulation and experimental data\n'
                        'example usage: python Get_Rates_Make_Trajectory.py '
                        '-PlotSim')
    parser.add_argument('-CreateTikzData',action='store_true')
    parser.add_argument('-PlotCoarse',action='store_true')

    args = parser.parse_args()

    MultiBodyCase = False

    if args.ExpData:

        if len(args.ExpData) > 1:
            MultiBodyCase = True

        if MultiBodyCase:

            Nbody = len(args.ExpData)
            # We create a list of classes here.. this is really cool.
            experimentalData = [Experimental_Dynamic_Data(args.ExpData[i]) for i in range(Nbody)] # Holy fuck this is awesome


            # Now we will do the typical options for single body, but for each object
            for i in range(Nbody):

                experimentalData[i].startup_MultiBody()
                experimentalData[i].Extract_Raw_Data(0,1,2,3,4,5,6)

                if args.mm2m:
                    experimentalData[i].X_loc = experimentalData[i].X_loc/1000
                    experimentalData[i].Y_loc = experimentalData[i].Y_loc/1000
                    experimentalData[i].Z_loc = experimentalData[i].Z_loc/1000

                if args.rad2deg:
                    experimentalData[i].Pitch = experimentalData[i].Pitch*180.0/np.pi
                    experimentalData[i].Yaw = experimentalData[i].Yaw*180.0/np.pi
                    experimentalData[i].Roll = experimentalData[i].Roll*180.0/np.pi

                if args.PlotExp:
                    plt.title("Body%i Raw Exp Data"%(i+1))
                    plt.plot(experimentalData[i].Time, experimentalData[i].X_loc, 'ro')
                    plt.plot(experimentalData[i].Time, experimentalData[i].Y_loc, 'bo')
                    plt.plot(experimentalData[i].Time, experimentalData[i].Z_loc, 'go')
                    plt.show()


                if args.PlotFits:

                    # Defining some standard/typical values
                    NFit = 500
                    XYZFitType = 'slinear'
                    PitchFitType = 'slinear'
                    YawFitType = 'slinear'

                    if args.PitchFitType:
                        PitchFitType = args.PitchFitType

                    if args.YawFitType:
                        YawFitType = args.YawFitType

                    if args.XYZFitType:
                        XYZFitType = args.XYZFitType

                    if args.NFit:
                        NFit = args.NFit

                    experimentalData[i].functional_fit(XYZFitType,PitchFitType,YawFitType,NFit)
                    experimentalData[i].Compute_X_VoT()

                    plt.title("Body%i Time vs. X"%(i+1))
                    plt.plot( experimentalData[i].Time,     experimentalData[i].X_loc, 'o', label='X-location' )
                    plt.plot( experimentalData[i].Time_Fit, experimentalData[i].X_loc_Fit, 'r-',  label='X_Fit-location' )
                    plt.legend()
                    plt.show()
                    plt.title("Body%i Time vs. X-V0T"%(i+1))
                    plt.plot( experimentalData[i].Time_Fit, experimentalData[i].X_VoT, 'r-', label='X-V0T' )
                    plt.legend()
                    plt.show()
                    plt.title("Body%i Time vs Y/Z"%(i+1))
                    plt.plot( experimentalData[i].Time,     experimentalData[i].Y_loc, 'o', label='Y-location' )
                    plt.plot( experimentalData[i].Time_Fit, experimentalData[i].Y_loc_Fit, 'r-', label='y_Fit-location' )
                    plt.plot( experimentalData[i].Time,     experimentalData[i].Z_loc, 'o', label='Z-location' )
                    plt.plot( experimentalData[i].Time_Fit, experimentalData[i].Z_loc_Fit, 'r-', label='Z_Fit-location' )
                    plt.legend()
                    plt.show()
                    plt.title("Body%i Time vs. Pitch"%(i+1))
                    plt.plot( experimentalData[i].Time,     experimentalData[i].Pitch, 'o', label='Pitch' )
                    plt.plot( experimentalData[i].Time_Fit, experimentalData[i].Pitch_Fit, 'r-', label='Pitch_Fit' )
                    plt.legend()
                    plt.show()
                    plt.title("Body%i Time vs. Yaw"%(i+1))
                    plt.plot( experimentalData[i].Time,     experimentalData[i].Yaw, 'o', label='Yaw' )
                    plt.plot( experimentalData[i].Time_Fit, experimentalData[i].Yaw_Fit, 'r-', label='Yaw_Fit' )
                    plt.legend()
                    plt.show()

                    experimentalData[i].estimate_rates()

            if args.PlotFits:
                plt.title("All Bodies XZ")
                for i in range(Nbody):
                    if i==0:
                        dx = experimentalData[i].X_locRaw[0]
                        dy = experimentalData[i].Y_locRaw[0]
                        dz = experimentalData[i].Z_locRaw[0]
                        plt.plot( experimentalData[i].X_loc_Fit, experimentalData[i].Z_loc_Fit, '-', label='X-Z-Fit Body%i'%(i+1) )
                    else:
                        dx = experimentalData[i].X_locRaw[0] - experimentalData[0].X_locRaw[0]
                        dy = experimentalData[i].Y_locRaw[0] - experimentalData[0].Y_locRaw[0]
                        dz = experimentalData[i].Z_locRaw[0] - experimentalData[0].Z_locRaw[0]
                        if args.mm2m:
                            dx = dx/1000
                            dy = dy/1000
                            dz = dz/1000
                        plt.plot( experimentalData[i].X_loc_Fit+dx, experimentalData[i].Z_loc_Fit+dz, '-', label='X-Z-Fit Body%i'%(i+1) )

                        # plt.plot( experimentalData[i].X_loc_Fit-experimentalData[0].X_loc_Fit, experimentalData[i].Z_loc_Fit, '-', label='X-Z-Fit Body%i'%(i+1) )

                plt.legend()
                # plt.xticks(np.arange(0.0,0.3, 0.01))
                # plt.yticks(np.arange(-0.5,0.0, 0.01))
                plt.grid()
                plt.show()

            if args.CreateTrajectory:
                
                # Defining some standard/typical values
                NFit = 500
                XYZFitType = 'slinear'
                PitchFitType = 'slinear'
                YawFitType = 'slinear'

                if args.PitchFitType:
                    PitchFitType = args.PitchFitType

                if args.YawFitType:
                    YawFitType = args.YawFitType

                if args.XYZFitType:
                    XYZFitType = args.XYZFitType

                if args.NFit:
                    NFit = args.NFit

                TrajFile = "Trajectory.dat"
                WriteHeader = True
                for j in range(Nbody):
                    print("\nFinding Fit data for %s"%experimentalData[j].dfile)
                    experimentalData[j].functional_fit(XYZFitType,PitchFitType,YawFitType,NFit)
                    experimentalData[j].estimate_rates()
                    experimentalData[j].Compute_X_VoT()


                # Create a single body trajectory file
                # experimentalData.Write_TrajectoryFile(TrajFile,WriteHeader,'Euler')

                # Create a multi body trajectory file but with frame discplacement
                print("Creating Trajectory File %s" % (TrajFile))
                for k in range(2):
                    # Doing the 'Type' thing now since this will be moved to another routine later...
                    if k == 0:
                        Type = 'Euler'
                    else:
                        Type = 'Lagrange'
                        TrajFile = TrajFile+"_Lagrange"
                    f1 = open(str(TrajFile),"wt")

                    if WriteHeader:
                        for j in range(Nbody):
                            if j==0:
                                f1.write('#Variables="time",')
                            f1.write('"xloc%i",'  % (j+1) )
                            f1.write('"yloc%i",'  % (j+1) )
                            f1.write('"zloc%i",'  % (j+1) )
                            f1.write('"theta%i",' % (j+1) )
                            f1.write('"psi%i",'   % (j+1) )
                            if j == Nbody - 1:
                                f1.write('"phi%i"\n' % (j+1) )
                            else:
                                f1.write('"phi%i",' % (j+1) )
                    for i in range(experimentalData[j].NFit):
                        for j in range(Nbody):
                            f1.write("%.9f      " %experimentalData[j].Time_Fit[i])
                            if Type == 'Lagrange':
                                f1.write("%.9f      " %experimentalData[j].X_VoT[i])
                            if Type == 'Euler':
                                f1.write("%.9f      " %experimentalData[j].X_loc_Fit[i])
                            f1.write("%.9f      " %experimentalData[j].Y_loc_Fit[i])
                            f1.write("%.9f      " %experimentalData[j].Z_loc_Fit[i])
                            f1.write("%.9f      " %experimentalData[j].Pitch_Fit[i])
                            f1.write("%.9f      " %experimentalData[j].Yaw_Fit[i])
                            if j == Nbody - 1:
                                f1.write("%.9f\n" %0.0) # Keeping this zero for now
                            else:
                                f1.write("%.9f      " %0.0) # Keeping this zero for now
                sys.exit("DEBUG3")

                f1.close()

                # Create a single body trajectory file but with frame discplacement
                experimentalData.Write_TrajectoryFile(TrajFile+"_Lagrange",WriteHeader,'Lagrange')

            sys.exit("DEBUG")

        else:
            ExpDataFile = args.ExpData[0]
            experimentalData = Experimental_Dynamic_Data(ExpDataFile)
            experimentalData.startup()
            # experimentalData.Extract_Raw_Data(0,3,1,2,4,5,6) # SIAD?
            experimentalData.Extract_Raw_Data(1,4,2,3,6,5,-1) # MPCV
            # experimentalData.Extract_Raw_Data(1,4,5,6,9,8,-1) # ADEPT
            # experimentalData.Extract_Raw_Data(0,1,2,3,5,4,-1) # MSL Aberdeen

            if args.mm2m:
                experimentalData.X_loc = experimentalData.X_loc/1000
                experimentalData.Y_loc = experimentalData.Y_loc/1000
                experimentalData.Z_loc = experimentalData.Z_loc/1000

            if args.rad2deg:
                experimentalData.Pitch = experimentalData.Pitch*180.0/np.pi
                experimentalData.Yaw = experimentalData.Yaw*180.0/np.pi
                experimentalData.Roll = experimentalData.Roll*180.0/np.pi


            if args.CreateTikzData:
                experimentalData.TikzData('exp_data.dat')


            if args.PlotExp:
                plt.plot(experimentalData.Time, experimentalData.Pitch, 'ro-',label="Pitch")
                plt.plot(experimentalData.Time, experimentalData.Yaw, 'bo-',label="Yaw")
                plt.legend()
                plt.show()

            if args.PlotFits:

                # Defining some standard/typical values
                NFit = 500
                XYZFitType = 'slinear'
                PitchFitType = 'slinear'
                YawFitType = 'slinear'

                if args.PitchFitType:
                    PitchFitType = args.PitchFitType

                if args.YawFitType:
                    YawFitType = args.YawFitType

                if args.XYZFitType:
                    XYZFitType = args.XYZFitType

                if args.NFit:
                    NFit = args.NFit

                experimentalData.functional_fit(XYZFitType,PitchFitType,YawFitType,NFit)

                plt.plot( experimentalData.Time, experimentalData.X_loc, 'o', label='X-location' )
                plt.plot( experimentalData.Time_Fit, experimentalData.X_loc_Fit, 'r-',  label='X_Fit-location' )
                plt.legend()
                plt.show()
                plt.plot( experimentalData.Time, experimentalData.Y_loc, 'o', label='Y-location' )
                plt.plot( experimentalData.Time_Fit, experimentalData.Y_loc_Fit, 'r-', label='y_Fit-location' )
                plt.plot( experimentalData.Time, experimentalData.Z_loc, 'o', label='Z-location' )
                plt.plot( experimentalData.Time_Fit, experimentalData.Z_loc_Fit, 'r-', label='Z_Fit-location' )
                plt.legend()
                plt.show()
                plt.plot( experimentalData.Time, experimentalData.Pitch, 'o', label='Pitch' )
                plt.plot( experimentalData.Time_Fit, experimentalData.Pitch_Fit, 'r-', label='Pitch_Fit' )
                plt.legend()
                plt.show()
                plt.plot( experimentalData.Time, experimentalData.Yaw, 'o', label='Yaw' )
                plt.plot( experimentalData.Time_Fit, experimentalData.Yaw_Fit, 'r-', label='Yaw_Fit' )
                plt.legend()
                plt.show()
                plt.plot( experimentalData.Time, experimentalData.Pitch, 'o', label='Pitch' )
                plt.plot( experimentalData.Time_Fit, experimentalData.Pitch_Fit, 'r-', label='Pitch_Fit' )
                plt.plot( experimentalData.Time, experimentalData.Yaw, 'o', label='Yaw' )
                plt.plot( experimentalData.Time_Fit, experimentalData.Yaw_Fit, 'r-', label='Yaw_Fit' )
                plt.legend()
                plt.show()

                experimentalData.estimate_rates()

                if args.CreateTrajectory:
                    TrajFile = "Trajectory.dat"
                    WriteHeader = True
                    # Create a single body trajectory file
                    experimentalData.Write_TrajectoryFile(TrajFile,WriteHeader,'Euler')
                    experimentalData.Compute_X_VoT()
                    # Create a single body trajectory file but with frame discplacement
                    experimentalData.Write_TrajectoryFile(TrajFile+"_Lagrange",WriteHeader,'Lagrange')
            
    if args.SimData:
        SimDataFile = args.SimData

        simulationData = Simulation_Dynamic_Data(SimDataFile)
        if args.RepairSimData:
            simulationData.readBadData()
            simulationData.readData()
        else:
            simulationData.readData()

        simulationData.ExtractData_SingleBody(0,1,2,3,16,17,18) # RigidBody.dat
        # simulationData.ExtractData_SingleBody(0,1,2,3,10,11,12) # dynamics_log.dat


        if args.PlotSim:
            plt.plot(simulationData.Time, simulationData.Pitch, 'r-') 
            plt.plot(simulationData.Time, simulationData.Yaw, 'b-')
            plt.show()

        if args.PlotCompare:
            if not args.ExpData:
                print("Need Experimental data!!")
                sys.exit("Exiting...")
            plt.plot(simulationData.Time, simulationData.Pitch, 'r-')
            plt.plot(simulationData.Time, simulationData.Yaw, 'b-')
            plt.plot(experimentalData.Time, experimentalData.Pitch, 'ro')
            plt.plot(experimentalData.Time, experimentalData.Yaw, 'bo')
            plt.show()

        if args.PlotCoarse:
            if not args.ExpData:
                print("Need Experimental data!!")
                sys.exit("Exiting...")
            simulationData.CoarsenData(int(simulationData.Ndata/500))
            plt.plot(simulationData.TimeC, simulationData.PitchC, 'r-')
            plt.plot(simulationData.TimeC, simulationData.YawC, 'b-')
            plt.plot(experimentalData.Time, experimentalData.Pitch, 'ro')
            plt.plot(experimentalData.Time, experimentalData.Yaw, 'bo')
            plt.show()

            experimentalData.Compute_X_VoT()
            simulationData.Compute_X_VoT()
            plt.title("X-V0T")
            plt.plot(experimentalData.Time_Fit, experimentalData.X_VoT, 'ro',label="Experiment")
            plt.plot(simulationData.TimeC, -simulationData.X_locC, 'b-',label='Simulation')
            plt.show()

        if args.CreateTikzData:
            simulationData.TikzData('sim_data.dat',int(simulationData.Ndata/500))

    print("Program complete!")
if __name__ == "__main__":
    main(sys.argv[1:])
