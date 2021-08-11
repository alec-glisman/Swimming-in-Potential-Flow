# External Dependencies
import json                        # Load and parse simulation metadata
import os                          # Access system file-tree
import sys                         # Modify system parameters
import matplotlib.pyplot as plt    # Plots (duh)
import numpy as np                 # Data structures
from optparse import OptionParser  # Get user input

# Internal Dependencies
sys.path.insert(0, os.getcwd() + '/python')
from plotStyling import PlotStyling  # noqa: E402


# NOTE: Parse user input
parser = OptionParser()
parser.add_option("--relPath", dest="user_relPath",
                  help="Confine bodies to only move in the x direction", metavar="PathString")
parser.add_option("--outputDir", dest="user_outputDir",
                  help="The ratio of the greater to lesser spring constants", metavar="PathString")


# Date that simulations of interest were run
def aggregate_plots(relPath, outputDir):

    # Parameters for function
    outputDir = relPath + outputDir + "/"
    metadata = []  # Lists to store loaded data
    epsOutput = True
    numPeriods = 1

    # REVIEW
    omega = 0.0

    # Load all data files
    try:
        # Loop through all subdirectories in the 'data' directory
        for root, dirs, files in os.walk(relPath, topdown=True):
            dirs.sort()  # Sort the directories

            # Loop through all files in given directory
            for file in files:

                # Only read data here from json files
                if ("json" in file):
                    with open(root + "/" + file) as j:
                        metadata.append(json.load(j))

        # Find \omega and U_0 from metadata output
        # omega   = float(metadata[0]["constraint"]["1D_kinematic"]["omega"])
        a = float(metadata[0]["scaling"]["length"])

    except:  # No files found
        print(f"No files found in relPath {relPath}")
        print(f"Failure to load data")
        raise IOError("FAILURE!!")

    # Loop through all simulations and grab the final CoM displacement (x-axis)
    CoM_disp_x = np.zeros(len(metadata))
    relDispEqbm = np.zeros(len(metadata))
    phaseShift = np.zeros(len(metadata))

    for i in range(len(metadata)):
        CoM_disp_x[i] = float(
            metadata[i]["KPIs"]["final_CoM_displacement_x"]) / (numPeriods * a)
        relDispEqbm[i] = float(metadata[i]["KPIs"]["rel_disp_eqbm"]) / (a)
        phaseShift[i] = float(metadata[i]["constraint"]
                              ["1D_kinematic"]["phase_angle"])

    # PLOT: net displacement of swimmer vs. distance between spheres
    CoM_Plot = PlotStyling(r"$|X_0 / a |$", r"$ \Delta Z / a $",
                           title=None, loglog=False,
                           outputDir=outputDir, figName="forced-potential-oscillation-CoM_x-disp", eps=epsOutput,
                           continuousColors=False)
    # Show numerical data points
    CoM_Plot.make_plot()
    CoM_Plot.add_scatter(relDispEqbm, CoM_disp_x, zorder=2, label="Simulation")
    # Add legend
    CoM_Plot.add_legend(title=None, loc='best',
                        bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    # Adjust ticks and tick labels
    CoM_Plot.ax.set_xlim([2, 6])
    CoM_Plot.set_major_minor_ticks(
        xMajorLoc=1, xMinorLoc=0.5, yMajorLoc=None, yMinorLoc=None)
    CoM_Plot.set_yaxis_scientific()
    CoM_Plot.save_plot()

    # PLOT: log-log of net displacement of swimmer vs. distance between spheres
    CoM_PlotLL = PlotStyling(r"$|X_0 / a |$", r"$| \Delta Z / a |$",
                             title=None, loglog=True,
                             outputDir=outputDir, figName="forced-potential-oscillation-CoM_x-disp-loglog", eps=epsOutput,
                             continuousColors=False)
    CoM_PlotLL.make_plot()
    CoM_PlotLL.add_scatter(relDispEqbm, np.abs(
        CoM_disp_x), zorder=2, label="Simulation")
    CoM_PlotLL.save_plot()

    # PLOT: net displacement of swimmer vs phase Shift
    phaseShift_Plot = PlotStyling(r"Phase Shift, $\delta$", r"$\Delta Z / a$",
                                  title=None, loglog=False,
                                  outputDir=outputDir, figName="forced-potential-oscillation-phaseShift", eps=epsOutput,
                                  continuousColors=False)
    phaseShift_Plot.make_plot()
    phaseShift_Plot.add_scatter(
        phaseShift, CoM_disp_x, zorder=2, label="Simulation")
    phaseShift_Plot.add_legend(
        title=None, loc='best', bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    phaseShift_Plot.set_yaxis_scientific()
    phaseShift_Plot.save_plot()

    # Output data to log file
    file = outputDir + "output.txt"
    with open(file, "w") as f:
        # Find scaling of CoM_dis_x vs. relDispEqbm the simulation data points
        with CoM_Plot.suppress_stdout():
            m, b = np.polyfit(np.log(relDispEqbm),
                              np.log(np.abs(CoM_disp_x)), 1)
        print(f"ln(CoM_disp_x) = {m} * ln(relDispEqbm) + {b}", file=f)
        print("\n", file=f)

        # Dump any other relevant data
        print(f"omega = {omega}",   file=f)
        print(f"a = {a}",           file=f)
        print("",                   file=f)
        print("phaseShift",      file=f)
        print(phaseShift,        file=f)
        print("relDispEqbm",     file=f)
        print(relDispEqbm,       file=f)
        print("CoM_disp_x",      file=f)
        print(CoM_disp_x,        file=f)
        print("\n",              file=f)


# For use when being called from C++
if __name__ == "__main__":
    options, remainder = parser.parse_args(sys.argv[1:])

    RelPath = str(options.user_relPath)
    OutputDir = str(options.user_outputDir)

    aggregate_plots(RelPath, OutputDir)
