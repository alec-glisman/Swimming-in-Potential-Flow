# SECTION: Depdencendies

# External Dependencies
import os                          # Access system file-tree
import sys                         # Modify system parameters
from math import isclose           # isclose()
import numpy as np                 # Data structures
from optparse import OptionParser  # Get user input
import matplotlib.ticker as mticker  # Scientific notation in labels
from matplotlib.ticker import FuncFormatter

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


# Function for D: derivative of mass matrix elements, dimensionless
def f_D(x):
    # In R/a coordinates
    numerator = 3.0 * (-68.0 + 93.0 * np.power(x, 3))
    denominator = x * np.square(17.0 - 18.0 * np.power(x, 3))
    return np.divide(numerator, denominator)


# Function for \Delta Z / a
def dZ_leadingOrder(phi, U0, omega, a, x):
    # In R/a coordinates
    return -(np.pi * np.sin(phi)) * np.square(U0 / (a * omega)) * f_D(x)


# Function to calculate relative error
def relErr(exact, approximate):
    return np.divide(np.abs(exact - approximate), np.abs(exact))

# Function to load data, perform analysis, and generate plots


def aggregate_plots(relative_path, output_dir):

    # SECTION: Parameters for function

    output_dir = relative_path + "/" + output_dir + "/"
    gsd_files = []
    epsOutput = True
    a = 1

    # REVIEW[epic=Future Features]: Move this adjustment into plotting style library
    # Correctly get scientific notation in text elements
    def scientific(x, pos):
        return '%0.2e' % x

    scientific_formatter = FuncFormatter(scientific)
    fmt = mticker.FuncFormatter(scientific_formatter)

# !SECTION (Parameters for function)


# SECTION: Load data

    try:
        # Loop through all subdirectories in the 'data' directory
        for root, dirs, files in os.walk(relative_path, topdown=True):
            dirs.sort()
            files.sort()

            for file in files:  # Loop through all files in given directory
                if (".gsd" in file):
                    with open(root + "/" + file) as g:
                        cur_gsd = GSDUtil(g.name, create_gsd=False)
                        gsd_files.append(cur_gsd)

        assert(len(gsd_files) > 0)

    except:  # No files found
        raise IOError(
            f"Failure to load data. No files found in relPath {relative_path}")

    # Loop through all simulations and grab the final CoM displacement (x-axis)
    CoM_disp_x = np.zeros(len(gsd_files))
    relDispEqbm = np.zeros(len(gsd_files))
    phaseShift = np.zeros(len(gsd_files))
    U0 = np.zeros(len(gsd_files))
    omega = np.zeros(len(gsd_files))
    epsilon = np.zeros(len(gsd_files))
    final_t = np.zeros(len(gsd_files))

    for i in range(len(gsd_files)):

        # Data from final frame
        gsd_files[i].snapshot = gsd_files[i].trajectory.read_frame(
            gsd_files[i].trajectory.file.nframes - 1)
        final_t[i] = float(gsd_files[i].snapshot.log['integrator/t'])
        # particle 1, x-coordinate
        CoM_disp_x[i] = float(gsd_files[i].snapshot.particles.position[1][0])

        # Data from header frame
        gsd_files[i].snapshot = gsd_files[i].trajectory.read_frame(0)
        relDispEqbm[i] = float(gsd_files[i].snapshot.log['swimmer/R_avg'])
        phaseShift[i] = float(gsd_files[i].snapshot.log['swimmer/phase_shift'])
        U0[i] = float(gsd_files[i].snapshot.log['swimmer/U0'])
        omega[i] = float(gsd_files[i].snapshot.log['swimmer/omega'])
        epsilon[i] = U0[i] / relDispEqbm[i] / omega[i]

        # Data from initial frame
        gsd_files[i].snapshot = gsd_files[i].trajectory.read_frame(1)
        CoM_disp_x[i] -= float(gsd_files[i].snapshot.particles.position[1][0])

# !SECTION (Load data)


# SECTION: Output data

    # Output data to log file
    file = output_dir + "output.txt"

    with open(file, "w") as f:

        Empty_Plot = PlotStyling(1, r"", r"",
                                 title=None, loglog=False,
                                 outputDir=output_dir, figName="", eps=epsOutput,
                                 continuousColors=False)

        # Find scaling of CoM_dis_x vs. relDispEqbm the simulation data points
        if (len(np.unique(relDispEqbm)) > 1):
            try:
                with Empty_Plot.suppress_stdout():
                    m, b = np.polyfit(np.log(np.abs(relDispEqbm)), np.log(
                        np.abs(CoM_disp_x), 1))
                    print(
                        f"ln(CoM_disp_x) = {m} * ln(relDispEqbm) + {b}",  file=f)
                    print("\n", file=f)
            except:
                pass

        print("CoM_disp_x:", file=f)
        print(CoM_disp_x,   file=f)
        print("", file=f)

        print("relDispEqbm:", file=f)
        print(relDispEqbm,    file=f)
        print("", file=f)

        print("final_t:", file=f)
        print(final_t,    file=f)

# !SECTION (Output data)


# SECTION: Analysis

    # Calculate leading order net motion over one period of articulation (for varying distance between spheres)
    xAnalyticalRng = np.array(np.linspace(
        2.0, 40.0, num=1000), dtype=np.float64)
    dZAnalyticalDist = dZ_leadingOrder(
        phaseShift[0], U0[0], omega[0], a, xAnalyticalRng)
    # Calculate leading order net motion over one period of articulation (for varying phase Shift)
    deltaAnalyticalRng = np.array(np.linspace(
        0, 2 * np.pi, num=1000), dtype=np.float64)
    dZAnalyticaldelt = dZ_leadingOrder(
        deltaAnalyticalRng, U0[0], omega[0], a, relDispEqbm[0])
    # Calculate leading order net motion (for varying epsilon)
    epsAnalyticalRng = np.array(np.linspace(
        0, np.max(epsilon), num=1000), dtype=np.float64)
    dZAnalyticaleps = dZ_leadingOrder(
        phaseShift[0], epsAnalyticalRng * omega[0] * relDispEqbm[0], omega[0], a, relDispEqbm[0])

# !SECTION (Analysis)


# SECTION: Plots

    if (len(np.unique(relDispEqbm)) > 1):
        # PLOT: net displacement of swimmer vs. distance between spheres
        numLines = 2
        CoM_Plot = PlotStyling(numLines,
                               r"$\mathrm{R}_0 / a $", r"$\Delta \mathrm{R}_{2} / a$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="collinear-swimmer-CoM_x-disp", eps=epsOutput,
                               continuousColors=False)
        # Show numerical data points
        CoM_Plot.make_plot()
        CoM_Plot.curve(np.abs(xAnalyticalRng),
                       dZAnalyticalDist, zorder=1, label="Leading Order")
        CoM_Plot.scatter(relDispEqbm, CoM_disp_x,
                         zorder=2, label="Simulation")
        # Add legend
        CoM_Plot.legend(title=r"$\epsilon \leq$" + "{}".format(fmt(np.max(epsilon))),
                        loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        # Adjust ticks and tick labels
        CoM_Plot.ax.set_xlim([1.9, 6])
        # CoM_Plot.set_major_minor_ticks(
        #     xMajorLoc=1, xMinorLoc=0.5, yMajorLoc=None, yMinorLoc=None)
        CoM_Plot.set_yaxis_scientific()
        CoM_Plot.save_plot()

        # PLOT: log-log of net displacement of swimmer vs. distance between spheres
        numLines = 2
        CoM_PlotLL = PlotStyling(numLines,
                                 r"$\mathrm{R}_0 / a $", r"$\Delta \mathrm{R}_{2} / a$",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-CoM_x-disp-loglog", eps=epsOutput,
                                 continuousColors=False)
        CoM_PlotLL.make_plot()
        CoM_PlotLL.curve(np.abs(xAnalyticalRng), np.abs(
            dZAnalyticalDist), zorder=1, label="Leading Order")
        CoM_PlotLL.scatter(relDispEqbm, np.abs(
            CoM_disp_x), zorder=2, label="Simulation")
        CoM_PlotLL.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        CoM_PlotLL.save_plot()

        # PLOT: Relative error of displacement with relDisp
        numLines = 1
        relDisErr = relErr(dZ_leadingOrder(phaseShift[0], U0[0], omega[0], a, relDispEqbm),
                           CoM_disp_x)
        CoMDispErr_Plot = PlotStyling(numLines,
                                      r"$\mathrm{R}_0 / a $", r"Relative Error",
                                      title=None, loglog=True,
                                      outputDir=output_dir, figName="collinear-swimmer-CoM_x-disp-error", eps=epsOutput,
                                      continuousColors=False)
        CoMDispErr_Plot.make_plot()
        CoMDispErr_Plot.scatter(
            relDispEqbm, relDisErr, zorder=1, label="Relative Error")
        CoMDispErr_Plot.save_plot()

    if (len(np.unique(phaseShift)) > 1):
        # PLOT: net displacement of swimmer vs phase Shift
        numLines = 2
        phaseShift_Plot = PlotStyling(numLines,
                                      r"Phase Shift, $\delta$", r"$\Delta \mathrm{R}_{2} / a$",
                                      title=None, loglog=False,
                                      outputDir=output_dir, figName="collinear-swimmer-phaseShift", eps=epsOutput,
                                      continuousColors=False)
        phaseShift_Plot.make_plot()
        phaseShift_Plot.curve(
            deltaAnalyticalRng, dZAnalyticaldelt, zorder=1, label="Leading Order")
        phaseShift_Plot.scatter(
            phaseShift, CoM_disp_x, zorder=2, label="Simulation")
        phaseShift_Plot.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        phaseShift_Plot.set_yaxis_scientific()
        phaseShift_Plot.save_plot()

        # PLOT: Relative error of displacement with delta
        relPhErr = relErr(dZ_leadingOrder(phaseShift, U0[0], omega[0], a, relDispEqbm[0]),
                          CoM_disp_x)
        numLines = 1
        phaseShiftErr_Plot = PlotStyling(numLines,
                                         r"Phase Shift, $\delta$", r"Relative Error",
                                         title=None, loglog=True,
                                         outputDir=output_dir, figName="collinear-swimmer-phaseShift-error", eps=epsOutput,
                                         continuousColors=False)
        phaseShiftErr_Plot.make_plot()
        phaseShiftErr_Plot.scatter(
            phaseShift, relPhErr, zorder=1, label="Relative Error")
        phaseShiftErr_Plot.save_plot()

    if (len(np.unique(epsilon)) > 1):
        # PLOT: net displacement of swimmer vs epsilon
        numLines = 2
        eps_Plot = PlotStyling(numLines,
                               r"$\epsilon = \frac{\mathrm{U}_0 / \omega}{\mathrm{R}_0}$", r"$\Delta \mathrm{R}_{2} / a$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="collinear-swimmer-eps-scaling-CoM_x-disp", eps=epsOutput,
                               continuousColors=False)
        eps_Plot.make_plot()
        eps_Plot.curve(epsAnalyticalRng, dZAnalyticaleps,
                       zorder=1, label="Leading Order")
        eps_Plot.scatter(epsilon, CoM_disp_x,
                         zorder=2, label="Simulation")
        eps_Plot.legend(loc='best', bbox_to_anchor=(0.05, 0.05, 0.5, 0.9))
        eps_Plot.save_plot()

        # PLOT: log-log net displacement of swimmer vs epsilon
        numLines = 2
        epsLL_Plot = PlotStyling(numLines,
                                 r"$\epsilon = \frac{\mathrm{U}_0 / \omega}{\mathrm{R}_0}$", r"$\Delta \mathrm{R}_{2} / a$",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-eps-scaling-CoM_x-disp-loglog", eps=epsOutput,
                                 continuousColors=False)
        epsLL_Plot.make_plot()
        epsLL_Plot.curve(epsAnalyticalRng, dZAnalyticaleps,
                         zorder=1, label="Leading Order")
        epsLL_Plot.scatter(epsilon, CoM_disp_x,
                           zorder=2, label="Simulation")
        epsLL_Plot.legend(
            loc='best', bbox_to_anchor=(0.05, 0.01, 0.5, 0.98))
        epsLL_Plot.save_plot()
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
