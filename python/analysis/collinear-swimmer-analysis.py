# SECTION: Depdencendies

# External Dependencies
import os                          # Access system file-tree
import sys                         # Modify system parameters
import numpy as np                 # Data structures
from optparse import OptionParser  # Get user input

# Internal Dependencies
sys.path.insert(0, os.getcwd() + '/python')
from plotStyling import PlotStyling  # noqa: E402
from GSDUtil import GSDUtil  # noqa: E402

# !SECTION (Dependencies)


# SECTION: User input options

parser = OptionParser()
parser.add_option("--relative-path", dest="u_relative_path",
                  help="Path to parent of simulation data",
                  metavar="string")
parser.add_option("--output-dir", dest="u_output_dir",
                  help="Directory to output plots (relative to relative-path)",
                  metavar="string")

# !SECTION: (User input options)


# Date that simulations of interest were run
def aggregate_plots(relative_path, output_dir):

    # SECTION: Parameters for function

    output_dir = relative_path + "/" + output_dir + "/"
    gsd_files = []
    epsOutput = True

# !SECTION (Parameters for function)


# SECTION: Load data

    try:
        # Loop through all subdirectories in the 'data' directory
        for root, dirs, files in os.walk(relative_path, topdown=True):
            dirs.sort()  # Sort the directories

            for file in files:  # Loop through all files in given directory
                if (".gsd" in file):
                    with open(root + "/" + file) as g:
                        gsd_files.append(GSDUtil(g, create_gsd=False))

            assert(len(gsd_files) > 0)

    except:  # No files found
        raise IOError(
            f"Failure to load data. No files found in relPath {relative_path}")

    # Loop through all simulations and grab the final CoM displacement (x-axis)
    CoM_disp_x = np.zeros(len(gsd_files))
    relDispEqbm = np.zeros(len(gsd_files))
    phaseShift = np.zeros(len(gsd_files))

    for i in range(len(gsd_files)):

        CoM_disp_x[i] = float(
            gsd_files[i].snapshot.particles.position[3 * 1]
        )
        relDispEqbm[i] = float(gsd_files[i].snapshot.log['swimmer/R_avg'])
        phaseShift[i] = float(gsd_files[i].snapshot.log['swimmer/phase_shift'])

# !SECTION (Load data)


# SECTION: Plots

    # PLOT: net displacement of swimmer vs. distance between spheres
    CoM_Plot = PlotStyling(r"$|X_0 / a |$", r"$ \Delta Z / a $",
                           title=None, loglog=False,
                           outputDir=output_dir, figName="collinear-swimmer-CoM_x-disp", eps=epsOutput,
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
                             outputDir=output_dir, figName="collinear-swimmer-CoM_x-disp-loglog", eps=epsOutput,
                             continuousColors=False)
    CoM_PlotLL.make_plot()
    CoM_PlotLL.add_scatter(relDispEqbm, np.abs(
        CoM_disp_x), zorder=2, label="Simulation")
    CoM_PlotLL.save_plot()

    # PLOT: net displacement of swimmer vs phase shift
    phaseShift_Plot = PlotStyling(r"Phase Shift, $\delta$", r"$\Delta Z / a$",
                                  title=None, loglog=False,
                                  outputDir=output_dir, figName="collinear-swimmer-phaseShift", eps=epsOutput,
                                  continuousColors=False)
    phaseShift_Plot.make_plot()
    phaseShift_Plot.add_scatter(
        phaseShift, CoM_disp_x, zorder=2, label="Simulation")
    phaseShift_Plot.add_legend(
        title=None, loc='best', bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    phaseShift_Plot.set_yaxis_scientific()
    phaseShift_Plot.save_plot()

# !SECTION (Plots)


# SECTION: For use when being called from command line
if __name__ == "__main__":

    # parse user input
    options, remainder = parser.parse_args(sys.argv[1:])
    relative_path = str(options.u_relative_path)
    output_dir = str(options.u_output_dir)

    # generate plots
    aggregate_plots(relative_path, output_dir)

# !SECTION (For use when being called from command line)
