# SECTION: Depdencendies

# External Dependencies
import os                          # Access system file-tree
import sys                         # Modify system parameters
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
            dirs.sort()  # Sort the directories

            for file in files:  # Loop through all files in given directory
                if (".gsd" in file):
                    with open(root + "/" + file) as g:
                        cur_gsd = GSDUtil(g.name, create_gsd=False)
                        gsd_files.append(cur_gsd)

        assert(len(gsd_files) == 1)
        gsd_file = gsd_files[0]

    except:  # No files found
        raise IOError(
            f"Failure to load data. Not exactly 1 file found in relPath {relative_path}")

    # Data from header frame
    gsd_file.snapshot = gsd_file.trajectory.read_frame(0)
    relDispEqbm = float(gsd_file.snapshot.log['swimmer/R_avg'])
    phaseShift = float(gsd_file.snapshot.log['swimmer/phase_shift'])
    U0 = float(gsd_file.snapshot.log['swimmer/U0'])
    omega = float(gsd_file.snapshot.log['swimmer/omega'])
    tau = float(gsd_file.snapshot.log['integrator/tau'])
    epsilon = U0 / relDispEqbm / omega

    # Initialize temporal data
    nframes = gsd_file.trajectory.file.nframes
    time = np.zeros((nframes - 1))
    positions = np.zeros((9, nframes - 1))
    velocities = np.zeros_like(positions)
    accelerations = np.zeros_like(positions)

    # Loop over all snapshots in GSD (skipping header with incorrect kinematics)
    for i in range(1, nframes-1):
        current_snapshot = gsd_file.trajectory.read_frame(i)

        time[i] = current_snapshot.log['integrator/t']
        positions[:, i-1] = current_snapshot.particles.position.flatten()
        velocities[:, i-1] = current_snapshot.particles.position.flatten()
        accelerations[:, i-1] = current_snapshot.particles.position.flatten()

# !SECTION (Load data)


# SECTION: Plots
    # PLOT: Relative oscillator displacement (x-axis)
    numLines = 2
    sprDis_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta x \, \omega / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="spring_disp-x", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    sprDis_Plot.make_plot()
    sprDis_Plot.curve(
        time, (positions[0, :] - positions[3, :] - relDispEqbm) * omega / U0, zorder=1, label=r"$1-2$")
    sprDis_Plot.curve(
        time, (positions[3, :] - positions[6, :] - relDispEqbm) * omega / U0, zorder=2, label=r"$2-3$")
    # Add legend
    sprDis_Plot.legend(loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
    sprDis_Plot.save_plot()

    # PLOT: spring velocity (x)
    numLines = 4
    sprVel_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta U / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="spring_vel-x", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    sprVel_Plot.make_plot()
    sprVel_Plot.curve(time, (velocities[0, :] - velocities[3, :]) / (
        U0), zorder=1, label=r"$2-3$ Simulation")
    sprVel_Plot.curve(time, (velocities[3, :] - velocities[6, :]) / (
        U0), zorder=2, label=r"$1-2$ Simulation")
    sprVel_Plot.curve(
        time, -np.cos(2 * np.pi * time), thin_curve=True, zorder=3, label=r"$1-2$ Leading Order")
    sprVel_Plot.curve(
        time, -np.cos(2 * np.pi * time + phaseShift), thin_curve=True, zorder=4, label=r"$2-3$ Leading Order")
    # Add legend
    sprVel_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    sprVel_Plot.save_plot()

    # PLOT: spring acceleration (x)
    numLines = 2
    sprAcc_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta \dot{U} / (U_0 \, \omega)$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="spring_acc-x", eps=epsOutput,
                              continuousColors=False)

    # Show numerical data points
    sprAcc_Plot.make_plot()
    sprAcc_Plot.curve(time, (accelerations[0, :] - accelerations[3, :]) / (
        U0 * omega), zorder=1, label=r"$1-2$")
    sprAcc_Plot.curve(time, (accelerations[3, :] - accelerations[6, :]) / (
        U0 * omega), zorder=2, label=r"$2-3$")
    # Add legend
    sprAcc_Plot.legend(loc='best', bbox_to_anchor=(0.01, 0.01, 0.98, 0.98))
    sprAcc_Plot.save_plot()

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
