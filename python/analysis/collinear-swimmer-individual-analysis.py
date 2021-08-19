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

    # Data from initial frame (not 0)
    gsd_file.snapshot = gsd_file.trajectory.read_frame(1)
    relDispEqbm = gsd_file.snapshot.log['swimmer/R_avg']
    phaseShift = gsd_file.snapshot.log['swimmer/phase_shift']
    U0 = gsd_file.snapshot.log['swimmer/U0']
    omega = gsd_file.snapshot.log['swimmer/omega']
    tau = gsd_file.snapshot.log['integrator/tau']
    epsilon = U0 / relDispEqbm / omega

    # Initialize temporal data
    nframes = gsd_file.trajectory.file.nframes
    time = np.zeros((nframes - 2), dtype=np.float64)
    positions = np.zeros((9, nframes - 2), dtype=np.float64)
    velocities = np.zeros_like(positions, dtype=np.float64)
    accelerations = np.zeros_like(positions, dtype=np.float64)

    # Loop over all snapshots in GSD (skipping header with incorrect kinematics)
    for i in range(1, nframes-1):
        current_snapshot = gsd_file.trajectory.read_frame(i)

        time[i-1] = current_snapshot.log['integrator/t']
        positions[:, i-1] = current_snapshot.log['particles/double_position'].flatten()
        velocities[:, i-1] = current_snapshot.log['particles/double_velocity'].flatten()
        accelerations[:, i -
                      1] = current_snapshot.log['particles/double_moment_inertia'].flatten()

# !SECTION (Load data)


# SECTION: Output data

    # Output data to log file
    file = output_dir + "output.txt"
    with open(file, "w") as f:
        print(f"R_avg = {relDispEqbm}", file=f)
        print(f"\Delta R_2 = {positions[3, -1] - positions[3, 0]}", file=f)

# !SECTION (Output data)


# SECTION: Plots
    # PLOT: Locater point position (x-axis)
    numLines = 1
    locPos_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta R_2$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="loc-pos-x", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    locPos_Plot.make_plot()
    locPos_Plot.curve(
        time, positions[3, :], zorder=1, label=r"$2$")
    locPos_Plot.save_plot()

    # PLOT: Relative oscillator displacement (x-axis)
    numLines = 4
    oscDis_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta x \, \omega / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-disp-x", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    oscDis_Plot.make_plot()
    oscDis_Plot.curve(
        time, (positions[0, :] - positions[3, :] - relDispEqbm) * omega / U0, zorder=1, label=r"$1-2$")
    oscDis_Plot.curve(
        time, (positions[3, :] - positions[6, :] - relDispEqbm) * omega / U0, zorder=2, label=r"$2-3$")
    oscDis_Plot.curve(
        time, np.sin(omega * tau * time), thin_curve=True, zorder=3, label=r"$1-2$ Constraint")
    oscDis_Plot.curve(
        time, -np.sin(omega * tau * time + phaseShift), thin_curve=True, zorder=4, label=r"$2-3$ Constraint")
    # Add legend
    oscDis_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscDis_Plot.save_plot()

    numLines = 2
    oscDis_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"Error $\Delta x \, \omega / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-disp-x-err", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    oscDis_Plot.make_plot()
    oscDis_Plot.curve(
        time, (positions[0, :] - positions[3, :] - relDispEqbm) * omega / U0 - np.sin(omega * tau * time), zorder=1, label=r"$1-2$")
    oscDis_Plot.curve(
        time, (positions[3, :] - positions[6, :] - relDispEqbm) * omega / U0 + np.sin(omega * tau * time + phaseShift), zorder=2, label=r"$2-3$")
    # Add legend
    oscDis_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscDis_Plot.save_plot()

    # PLOT: Relative oscillator velocity (x-axis)
    numLines = 4
    oscVel_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta U / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-vel-x", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    oscVel_Plot.make_plot()
    oscVel_Plot.curve(time, (velocities[0, :] - velocities[3, :]) / (
        U0), zorder=1, label=r"$1-2$ Simulation")
    oscVel_Plot.curve(time, (velocities[6, :] - velocities[3, :]) / (
        U0), zorder=2, label=r"$3-2$ Simulation")
    oscVel_Plot.curve(
        time, np.cos(omega * tau * time), thin_curve=True, zorder=3, label=r"$1-2$ Constraint")
    oscVel_Plot.curve(
        time, np.cos(omega * tau * time + phaseShift), thin_curve=True, zorder=4, label=r"$3-2$ Constraint")
    # Add legend
    oscVel_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscVel_Plot.save_plot()

    numLines = 2
    oscVel_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta U / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-vel-x-err", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    oscVel_Plot.make_plot()
    oscVel_Plot.curve(time, (velocities[0, :] - velocities[3, :]) / (
        U0) - np.cos(omega * tau * time), zorder=1, label=r"$1-2$ Simulation")
    oscVel_Plot.curve(time, (velocities[6, :] - velocities[3, :]) / (
        U0) - np.cos(omega * tau * time + phaseShift), zorder=2, label=r"$3-2$ Simulation")
    # Add legend
    oscVel_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscVel_Plot.save_plot()

    # PLOT: Relative oscillator acceleration (x-axis)
    numLines = 4
    oscAcc_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta \dot{U} / (U_0 \, \omega)$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-acc-x", eps=epsOutput,
                              continuousColors=False)

    # Show numerical data points
    oscAcc_Plot.make_plot()
    oscAcc_Plot.curve(time, (accelerations[0, :] - accelerations[3, :]) / (
        U0 * omega), zorder=1, label=r"$1-2$ Simulation")
    oscAcc_Plot.curve(time, (accelerations[6, :] - accelerations[3, :]) / (
        U0 * omega), zorder=2, label=r"$3-2$ Simulation")
    oscAcc_Plot.curve(
        time, -np.sin(omega * tau * time), thin_curve=True, zorder=3, label=r"$1-2$ Constraint")
    oscAcc_Plot.curve(
        time, -np.sin(omega * tau * time + phaseShift), thin_curve=True, zorder=4, label=r"$3-2$ Constraint")
    # Add legend
    oscAcc_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscAcc_Plot.save_plot()

    numLines = 2
    oscAcc_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta \dot{U} / (U_0 \, \omega)$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-acc-x-err", eps=epsOutput,
                              continuousColors=False)

    # Show numerical data points
    oscAcc_Plot.make_plot()
    oscAcc_Plot.curve(time, (accelerations[0, :] - accelerations[3, :]) / (
        U0 * omega) + np.sin(omega * tau * time), zorder=1, label=r"$1-2$ Simulation")
    oscAcc_Plot.curve(time, (accelerations[6, :] - accelerations[3, :]) / (
        U0 * omega) + np.sin(omega * tau * time + phaseShift), zorder=2, label=r"$3-2$ Simulation")
    # Add legend
    oscAcc_Plot.legend(
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscAcc_Plot.save_plot()

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
