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
    CoM_disp_comp = np.zeros((len(gsd_files), 3), dtype=np.double)
    CoM_disp = np.zeros(len(gsd_files), dtype=np.double)
    R_avg = np.zeros_like(CoM_disp, dtype=np.double)
    Z_height = np.zeros_like(CoM_disp, dtype=np.double)
    phaseShift = np.zeros_like(CoM_disp, dtype=np.double)
    U0 = np.zeros_like(CoM_disp, dtype=np.double)
    omega = np.zeros_like(CoM_disp, dtype=np.double)
    epsilon = np.zeros_like(CoM_disp, dtype=np.double)
    final_t = np.zeros_like(CoM_disp, dtype=np.double)
    dt = np.zeros_like(CoM_disp, dtype=np.double)

    for i in range(len(gsd_files)):

        # Data from final frame
        gsd_files[i].snapshot = gsd_files[i].trajectory.read_frame(
            gsd_files[i].trajectory.file.nframes - 1)
        final_t[i] = gsd_files[i].snapshot.log['integrator/t']
        dt[i] = gsd_files[i].snapshot.log['integrator/dt']
        CoM_disp_comp[i, :] = gsd_files[i].snapshot.log['particles/double_position'][1]

        # Data from initial frame (not 0)
        gsd_files[i].snapshot = gsd_files[i].trajectory.read_frame(1)
        R_avg[i] = gsd_files[i].snapshot.log['swimmer/R_avg']
        Z_height[i] = gsd_files[i].snapshot.log['swimmer/Z_height']
        phaseShift[i] = gsd_files[i].snapshot.log['swimmer/phase_shift']
        U0[i] = gsd_files[i].snapshot.log['swimmer/U0']
        omega[i] = gsd_files[i].snapshot.log['swimmer/omega']
        epsilon[i] = U0[i] / R_avg[i] / omega[i]
        CoM_disp_comp[i, :] -= gsd_files[i].snapshot.log['particles/double_position'][1]

        CoM_disp[i] = np.linalg.norm(CoM_disp_comp[i, :])

# !SECTION (Load data)


# SECTION: Analysis

# !SECTION (Analysis)


# SECTION: Plots

    if (len(np.unique(R_avg)) > 1):

        idx = np.argsort(R_avg)
        CoM_disp_srt = CoM_disp[idx]
        R_avg_srt = R_avg[idx]

        # PLOT: net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_Plot = PlotStyling(numLines,
                               r"$\mathrm{R}_0 / a $", r"$\Delta \mathrm{R}_{2} / a$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="collinear-swimmer-wall-CoM_x-disp", eps=epsOutput,
                               continuousColors=False)
        # Show numerical data points
        CoM_Plot.make_plot()
        CoM_Plot.scatter_dashed(R_avg_srt[R_avg_srt <= 6.0], CoM_disp_srt[R_avg_srt <= 6.0],
                                zorder=1, label="Simulation")
        # Add legend
        CoM_Plot.legend(title=r"$\epsilon \leq$" + "{}".format(fmt(np.max(epsilon))),
                        loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        # Adjust ticks and tick labels
        # CoM_Plot.set_major_minor_ticks(
        #     xMajorLoc=1, xMinorLoc=0.5, yMajorLoc=None, yMinorLoc=None)
        CoM_Plot.set_yaxis_scientific()
        CoM_Plot.save_plot()

        # PLOT: log-log of net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_PlotLL = PlotStyling(numLines,
                                 r"$\mathrm{R}_0 / a $", r"$\Delta \mathrm{R}_{2} / a$",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-wall-CoM-disp-loglog", eps=epsOutput,
                                 continuousColors=False)
        CoM_PlotLL.make_plot()
        CoM_PlotLL.scatter_dashed(R_avg_srt, np.abs(
            CoM_disp_srt), zorder=2, label="Simulation")
        CoM_PlotLL.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        CoM_PlotLL.save_plot()

    if (len(np.unique(dt)) > 1):

        idx = np.argsort(dt)
        CoM_disp_srt = CoM_disp[idx]
        dt_srt = dt[idx]

        # PLOT: net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_Plot = PlotStyling(numLines,
                               r"$\Delta t / \tau$", r"$\Delta \mathrm{R}_{2} / a$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="collinear-swimmer-wall-CoM-disp-dtVary", eps=epsOutput,
                               continuousColors=False)
        # Show numerical data points
        CoM_Plot.make_plot()
        CoM_Plot.scatter_dashed(dt_srt, CoM_disp_srt,
                                zorder=1, label="Simulation")
        # Add legend
        CoM_Plot.legend(title=r"$\epsilon \leq$" + "{}".format(fmt(np.max(epsilon))),
                        loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        # Adjust ticks and tick labels
        # CoM_Plot.set_major_minor_ticks(
        #     xMajorLoc=1, xMinorLoc=0.5, yMajorLoc=None, yMinorLoc=None)
        CoM_Plot.set_yaxis_scientific()
        CoM_Plot.save_plot()

        # PLOT: log-log of net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_PlotLL = PlotStyling(numLines,
                                 r"$\Delta t / \tau$", r"$\Delta \mathrm{R}_{2} / a$",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-wall-CoM-disp-loglog-dtVary", eps=epsOutput,
                                 continuousColors=False)
        CoM_PlotLL.make_plot()
        CoM_PlotLL.scatter_dashed(dt_srt, np.abs(
            CoM_disp_srt), zorder=2, label="Simulation")
        CoM_PlotLL.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        CoM_PlotLL.save_plot()

        # PLOT: log-log of net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_PlotLL = PlotStyling(numLines,
                                 r"$\Delta t / \tau$", r"$\Delta \mathrm{R}_{2} / a - ($ val @ $\mathrm{min}(dt)$)",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-wall-CoM-disp-loglog-dtVary-diffFromLast", eps=epsOutput,
                                 continuousColors=False)
        CoM_PlotLL.make_plot()
        CoM_PlotLL.scatter_dashed(dt_srt, np.abs(
            CoM_disp_srt - CoM_disp_srt[0]), zorder=2, label="Simulation")
        CoM_PlotLL.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        CoM_PlotLL.save_plot()

    if (len(np.unique(Z_height)) > 1):

        idx = np.argsort(Z_height)
        CoM_disp_srt = CoM_disp[idx]
        Z_height_srt = Z_height[idx]

        # PLOT: net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_Plot = PlotStyling(numLines,
                               r"$\mathrm{Z}_0 / a $", r"$\Delta \mathrm{R}_{2} / a$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="collinear-swimmer-wall-CoM-disp-height", eps=epsOutput,
                               continuousColors=False)
        # Show numerical data points
        CoM_Plot.make_plot()
        CoM_Plot.scatter_dashed(Z_height_srt, CoM_disp_srt,
                                zorder=1, label="Simulation")
        # Add legend
        CoM_Plot.legend(title=r"$\epsilon \leq$" + "{}".format(fmt(np.max(epsilon))),
                        loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        # Adjust ticks and tick labels
        # CoM_Plot.set_major_minor_ticks(
        #     xMajorLoc=1, xMinorLoc=0.5, yMajorLoc=None, yMinorLoc=None)
        CoM_Plot.set_yaxis_scientific()
        CoM_Plot.save_plot()

        # PLOT: log-log of net displacement of swimmer vs. distance between spheres
        numLines = 1
        CoM_PlotLL = PlotStyling(numLines,
                                 r"$\mathrm{Z}_0 / a $", r"$\Delta \mathrm{R}_{2} / a$",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-wall-CoM-disp-height-loglog", eps=epsOutput,
                                 continuousColors=False)
        CoM_PlotLL.make_plot()
        CoM_PlotLL.scatter_dashed(Z_height_srt, np.abs(
            CoM_disp_srt), zorder=2, label="Simulation")
        CoM_PlotLL.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        CoM_PlotLL.save_plot()

    if (len(np.unique(phaseShift)) > 1):

        idx = np.argsort(phaseShift)
        CoM_disp_srt = CoM_disp[idx]
        phaseShift_srt = phaseShift[idx]

        # PLOT: net displacement of swimmer vs phase Shift
        numLines = 1
        phaseShift_Plot = PlotStyling(numLines,
                                      r"Phase Shift, $\delta$", r"$\Delta \mathrm{R}_{2} / a$",
                                      title=None, loglog=False,
                                      outputDir=output_dir, figName="collinear-swimmer-wall-phaseShift", eps=epsOutput,
                                      continuousColors=False)
        phaseShift_Plot.make_plot()
        phaseShift_Plot.scatter_dashed(
            phaseShift_srt, CoM_disp_srt, zorder=2, label="Simulation")
        phaseShift_Plot.legend(title=r"$\epsilon \leq$" + "{}".format(
            fmt(np.max(epsilon))), loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
        phaseShift_Plot.set_yaxis_scientific()
        phaseShift_Plot.save_plot()

    if (len(np.unique(epsilon)) > 1):
        # PLOT: net displacement of swimmer vs epsilon
        numLines = 1
        eps_Plot = PlotStyling(numLines,
                               r"$\epsilon = \frac{\mathrm{U}_0 / \omega}{\mathrm{R}_0}$", r"$\Delta \mathrm{R}_{2} / a$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="collinear-swimmer-wall-eps-scaling-CoM-disp", eps=epsOutput,
                               continuousColors=False)
        eps_Plot.make_plot()
        eps_Plot.scatter(epsilon, CoM_disp,
                         zorder=2, label="Simulation")
        eps_Plot.legend(loc='best', bbox_to_anchor=(0.05, 0.05, 0.5, 0.9))
        eps_Plot.save_plot()

        # PLOT: log-log net displacement of swimmer vs epsilon
        numLines = 1
        epsLL_Plot = PlotStyling(numLines,
                                 r"$\epsilon = \frac{\mathrm{U}_0 / \omega}{\mathrm{R}_0}$", r"$\Delta \mathrm{R}_{2} / a$",
                                 title=None, loglog=True,
                                 outputDir=output_dir, figName="collinear-swimmer-wall-eps-scaling-CoM-disp-loglog", eps=epsOutput,
                                 continuousColors=False)
        epsLL_Plot.make_plot()
        epsLL_Plot.scatter(epsilon, CoM_disp,
                           zorder=2, label="Simulation")
        epsLL_Plot.legend(
            loc='best', bbox_to_anchor=(0.05, 0.01, 0.5, 0.98))
        epsLL_Plot.save_plot()
# !SECTION (Plots)


# SECTION: Output data

    # Output data to log file
    file = output_dir + "output.txt"

    with open(file, "w") as f:

        Empty_Plot = PlotStyling(1, r"", r"",
                                 title=None, loglog=False,
                                 outputDir=output_dir, figName="", eps=epsOutput,
                                 continuousColors=False)

        # scaling of CoM_disp vs. R_avg the simulation data points
        if (len(np.unique(R_avg)) > 1):
            try:
                with Empty_Plot.suppress_stdout():
                    idx = np.argsort(R_avg)
                    CoM_disp_srt = CoM_disp[idx]
                    R_avg_srt = R_avg[idx]

                    m, b = np.polyfit(np.log(np.abs(R_avg_srt)),
                                      np.log(np.abs(CoM_disp_srt)), 1)

                print(f"ln(CoM_disp) = {m} * ln(R_avg) + {b}",  file=f)
                print("\n", file=f)
            except:
                pass

        # scaling of CoM_disp vs. Z_height the simulation data points
        if (len(np.unique(R_avg)) > 1):
            try:
                with Empty_Plot.suppress_stdout():
                    idx = np.argsort(Z_height)
                    CoM_disp_srt = CoM_disp[idx]
                    Z_height_srt = Z_height[idx]

                    m, b = np.polyfit(np.log(np.abs(Z_height)),
                                      np.log(np.abs(CoM_disp_srt)), 1)

                print(f"ln(CoM_disp) = {m} * ln(Z_height) + {b}",  file=f)
                print("\n", file=f)
            except:
                pass

        print("CoM_disp_comp:", file=f)
        print(CoM_disp_comp,   file=f)
        print("", file=f)

        print("CoM_disp:", file=f)
        print(CoM_disp,   file=f)
        print("", file=f)

        print("R_avg:", file=f)
        print(R_avg,    file=f)
        print("", file=f)

        print("Z_height:", file=f)
        print(Z_height,    file=f)
        print("", file=f)

        print("final_t:", file=f)
        print(final_t,    file=f)
        print("", file=f)

        print("dt:", file=f)
        print(dt, file=f)

# !SECTION (Output data)


# SECTION: For use when being called from command line
if __name__ == "__main__":

    # parse user input
    options, remainder = parser.parse_args(sys.argv[1:])
    relative_path = str(options.u_relative_path)
    output_dir = str(options.u_output_dir)

    # generate plots
    aggregate_plots(relative_path, output_dir)

# !SECTION (For use when being called from command line)
