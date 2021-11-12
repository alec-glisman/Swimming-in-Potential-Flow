"""Script to analyze a single GSD file.

__author__ = "Alec Glisman"

Example:
    python3 python/analysis/collinear-swimmer-wall-individual-analysis.py --relative-path=temp/output --output-dir=temp/output/figures
"""

# SECTION: Depdencendies

# External Dependencies
import os                          # Access system file-tree
import sys
# Modify system parameters
from matplotlib.pyplot import axis
import numpy as np                 # Data structures
from optparse import OptionParser  # Get user input
import matplotlib.ticker as mticker  # Scientific notation in labels
from matplotlib.ticker import FuncFormatter

# Internal Dependencies
sys.path.insert(0, os.getcwd() + '/python')
from plotStyling import PlotStyling  # noqa: E402
from GSDUtilPy import GSDUtilPy  # noqa: E402

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
    """Script to analyze a single GSD file.

    Args:
        relative_path (str): path to GSD files to load
        output_dir (str): path to output of numerical analysis

    Example:
        python3 python/analysis/collinear-swimmer-wall-individual-analysis.py --relative-path=temp/output --output-dir=temp/output/figures
    """

    # SECTION: Parameters for function
    output_dir = relative_path + "/" + output_dir + "/"
    gsd_files = []
    epsOutput = False
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
                        cur_gsd = GSDUtilPy(g.name, create_gsd=False)
                        gsd_files.append(cur_gsd)

        assert(len(gsd_files) == 1)
        gsd_file = gsd_files[0]

    except:  # No files found
        raise IOError(
            f"Failure to load data. Not exactly 1 file found in relPath {relative_path}")

    # Data from initial frame
    gsd_file.snapshot = gsd_file.trajectory.read_frame(0)
    N_particles = gsd_file.snapshot.particles.N
    nframes = gsd_file.trajectory.file.nframes
    R_avg = gsd_file.snapshot.log['swimmer/R_avg']
    Z_height = gsd_file.snapshot.log['swimmer/Z_height']
    phaseShift = gsd_file.snapshot.log['swimmer/phase_shift']
    U0 = gsd_file.snapshot.log['swimmer/U0']
    omega = gsd_file.snapshot.log['swimmer/omega']
    tau = gsd_file.snapshot.log['integrator/tau']
    epsilon = abs(U0 / R_avg / omega)

    # Initialize temporal data
    time = np.zeros((nframes - 2), dtype=np.float64)

    positions = np.zeros((N_particles, 3, nframes - 2), dtype=np.float64)
    velocities = np.zeros_like(positions, dtype=np.float64)
    accelerations = np.zeros_like(positions, dtype=np.float64)

    E_locater = np.zeros_like(time, dtype=np.float64)
    E_locater_internal = np.zeros_like(time, dtype=np.float64)
    E_internal = np.zeros_like(time, dtype=np.float64)
    E_simple = np.zeros_like(time, dtype=np.float64)

    # Loop over all snapshots in GSD (skipping header with incorrect kinematics)
    for i in range(1, nframes-1):
        current_snapshot = gsd_file.trajectory.read_frame(i)

        time[i-1] = current_snapshot.log['integrator/t']

        positions[:, :, i-1] = current_snapshot.log['particles/double_position']
        velocities[:, :, i-1] = current_snapshot.log['particles/double_velocity']
        accelerations[:, :, i -
                      1] = current_snapshot.log['particles/double_moment_inertia']

        E_locater[i-1] = current_snapshot.log['hydrodynamics/E_locater']
        E_locater_internal[i -
                           1] = current_snapshot.log['hydrodynamics/E_locater_internal']
        E_internal[i-1] = current_snapshot.log['hydrodynamics/E_internal']
        E_simple[i-1] = current_snapshot.log['hydrodynamics/E_simple']

# !SECTION (Load data)


# SECTION: Analysis

    # relative separation between particle pairs
    R_loc = positions[0, :, :]  # (spatial_dimension, frame_number)
    R_loc_init = positions[0, :, 0]
    DR_loc = R_loc - R_loc_init[:, np.newaxis]
    R_12 = positions[1, :, :] - R_loc
    R_32 = positions[2, :, :] - R_loc
    # distance between particle pairs
    Dr_Loc = np.linalg.norm(DR_loc, axis=0)
    r_12 = np.linalg.norm(R_12, axis=0)
    r_32 = np.linalg.norm(R_32, axis=0)

    # swimmer orientation
    q_unnormalized = positions[1, :, :] - positions[2, :, :]
    q_norms = np.linalg.norm(q_unnormalized, axis=0)
    q = q_unnormalized
    for i in range(q.shape[1]):
        q[:, i] /= q_norms[i]

    # swimmer angular rotation from initial configuration
    cos_theta = np.dot(q.T, q[:, 0])
    theta = np.arccos(cos_theta)
    theta_dot = np.gradient(theta, time)
    theta_ddot = np.gradient(theta_dot, time)

    # relative velocities between particle pairs
    U_loc = velocities[0, :, :]
    U_1 = velocities[1, :, :]
    U_3 = velocities[2, :, :]
    # swimmer velocity constraints
    u_12_sim = (q.T @ (U_1 - U_loc)).diagonal()
    u_32_sim = (q.T @ (U_3 - U_loc)).diagonal()

    # relative accelerations between particle pairs
    A_loc = accelerations[0, :, :]
    A_1 = accelerations[1, :, :]
    A_3 = accelerations[2, :, :]
    # swimmer acceleration constraints
    a_12_sim = (q.T @ (A_1 - A_loc)).diagonal()
    a_32_sim = (q.T @ (A_3 - A_loc)).diagonal()

    # characteristic scales
    char_time = 1.0 / omega
    char_vel = U0
    char_len = char_vel * char_time
    char_acc = char_vel / char_time

    # kinematic constraints
    r_12_con_mag = char_len * np.sin(omega * tau * time)
    r_32_con_mag = char_len * np.sin(omega * tau * time + phaseShift)
    u_12_con_mag = char_vel * np.cos(omega * tau * time)
    u_32_con_mag = char_vel * np.cos(omega * tau * time + phaseShift)
    a_12_con_mag = - char_acc * np.sin(omega * tau * time)
    a_32_con_mag = - char_acc * np.sin(omega * tau * time + phaseShift)

    # Energetics
    E_total = E_locater + E_locater_internal + E_internal
    E0 = E_total[0]

# !SECTION (Analysis)


# SECTION: Plots
    # PLOT: Energy dynamics
    numLines = 5
    energy_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$E / E_{0}$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="energy", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    energy_Plot.make_plot()
    energy_Plot.curve(time, E_internal / E0, zorder=0,
                      label=r"$E_{\mathrm{int}}$")
    energy_Plot.curve(time, E_locater_internal / E0, zorder=1,
                      label=r"$E_{\mathrm{loc-int}}$")
    energy_Plot.curve(time, E_locater / E0, zorder=2,
                      label=r"$E_{\mathrm{loc}}$")
    energy_Plot.curve(time, E_simple / E0, zorder=3,
                      label=r"$E_{\mathrm{s}}$")
    energy_Plot.curve(time, E_total / E0, zorder=4, label=r"$E$")
    energy_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=1, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    energy_Plot.save_plot()

    # PLOT: Angular displacement
    numLines = 1
    angDisp_Plot = PlotStyling(numLines,
                               r"$t/\tau$", r"$\Delta \theta / (2 \pi)$",
                               title=None, loglog=False,
                               outputDir=output_dir, figName="angular-displacement", eps=epsOutput,
                               continuousColors=False)
    # Show numerical data points
    angDisp_Plot.make_plot()
    angDisp_Plot.curve(
        time, theta / (2.0 * np.pi), zorder=1, label=r"Collinear Swimmer")
    angDisp_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=1, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    angDisp_Plot.save_plot()

    # PLOT: Angular velocity
    numLines = 1
    angVel_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta \dot{\theta} / (2 \pi)$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="angular-velocity", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    angVel_Plot.make_plot()
    angVel_Plot.curve(
        time, theta_dot / (2.0 * np.pi), zorder=1, label=r"Collinear Swimmer")
    angVel_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=1, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    angVel_Plot.save_plot()

    # PLOT: Angular acceleration
    numLines = 1
    angAcc_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta \ddot{\theta} / (2 \pi)$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="angular-acceleration", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    angAcc_Plot.make_plot()
    angAcc_Plot.curve(
        time, theta_ddot / (2.0 * np.pi), zorder=1, label=r"Collinear Swimmer")
    angAcc_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=1, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    angAcc_Plot.save_plot()

    # PLOT: Locater point position
    numLines = 1
    locPos_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$\Delta \mathbf{R}_2 \cdot \mathbf{q}$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="loc-pos", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    locPos_Plot.make_plot()
    locPos_Plot.curve(
        time, Dr_Loc, zorder=1, label=r"$2$")
    locPos_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=1, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    locPos_Plot.save_plot()

    # PLOT: Relative oscillator displacement (x-axis)
    numLines = 4
    oscDis_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$(\mathbf{r} \cdot \mathbf{q}) \omega / U_0$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-disp", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    oscDis_Plot.make_plot()
    oscDis_Plot.curve(
        time, (r_12 - R_avg) / char_len, zorder=1, label=r"$1-2$ Simulation")
    oscDis_Plot.curve(
        time, (r_32 - R_avg) / char_len, zorder=2, label=r"$3-2$ Simulation")
    oscDis_Plot.curve(
        time, r_12_con_mag / char_len, thin_curve=True, zorder=3, label=r"$1-2$ Constraint")
    oscDis_Plot.curve(
        time, -r_32_con_mag / char_len, thin_curve=True, zorder=4, label=r"$3-2$ Constraint")
    # Add legend
    oscDis_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscDis_Plot.save_plot()

    numLines = 2
    oscDisErr_Plot = PlotStyling(numLines,
                                 r"$t/\tau$", r"Error $(\mathbf{r} \cdot \mathbf{q}) \omega / U_0$",
                                 title=None, loglog=False,
                                 outputDir=output_dir, figName="osc-disp-err", eps=epsOutput,
                                 continuousColors=False)
    # Show numerical data points
    oscDisErr_Plot.make_plot()
    oscDisErr_Plot.curve(
        time, (r_12 - R_avg) / char_len - (r_12_con_mag) / char_len, zorder=1, label=r"$1-2$")
    oscDisErr_Plot.curve(
        time, (r_32 - R_avg) / char_len + (r_32_con_mag) / char_len, zorder=2, label=r"$2-3$")
    oscDisErr_Plot.set_yaxis_scientific()
    # Add legend
    oscDisErr_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=2, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    oscDisErr_Plot.save_plot()

    # PLOT: Relative oscillator velocity (x-axis)
    numLines = 4
    oscVel_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$(\Delta \mathbf{U} / U_0) \cdot \mathbf{q}$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-vel", eps=epsOutput,
                              continuousColors=False)
    # Show numerical data points
    oscVel_Plot.make_plot()
    oscVel_Plot.curve(time,  u_12_sim / char_vel, zorder=1,
                      label=r"$1-2$ Simulation")
    oscVel_Plot.curve(time, u_32_sim / char_vel, zorder=2,
                      label=r"$3-2$ Simulation")
    oscVel_Plot.curve(
        time, u_12_con_mag / char_vel, thin_curve=True, zorder=3, label=r"$1-2$ Constraint")
    oscVel_Plot.curve(
        time, u_32_con_mag / char_vel, thin_curve=True, zorder=4, label=r"$3-2$ Constraint")
    # Add legend
    oscVel_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscVel_Plot.save_plot()

    numLines = 2
    oscVelErr_Plot = PlotStyling(numLines,
                                 r"$t/\tau$", r"Error $(\Delta \mathbf{U} / U_0) \cdot \mathbf{q}$",
                                 title=None, loglog=False,
                                 outputDir=output_dir, figName="osc-vel-err", eps=epsOutput,
                                 continuousColors=False)
    # Show numerical data points
    oscVelErr_Plot.make_plot()
    oscVelErr_Plot.curve(time, u_12_sim / char_vel -
                         u_12_con_mag / char_vel, zorder=1, label=r"$1-2$")
    oscVelErr_Plot.curve(time, u_32_sim / char_vel -
                         u_32_con_mag / char_vel, zorder=2, label=r"$3-2$")
    oscVelErr_Plot.set_yaxis_scientific()
    # Add legend
    oscVelErr_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=2, bbox_to_anchor=(0.05, 0.05, 0.9, 0.9))
    oscVelErr_Plot.save_plot()

    # PLOT: Relative oscillator acceleration (x-axis)
    numLines = 4
    oscAcc_Plot = PlotStyling(numLines,
                              r"$t/\tau$", r"$(\Delta \dot{\mathbf{U}} \cdot \mathbf{q}) / (U_0 \, \omega)$",
                              title=None, loglog=False,
                              outputDir=output_dir, figName="osc-acc", eps=epsOutput,
                              continuousColors=False)

    # Show numerical data points
    oscAcc_Plot.make_plot()
    oscAcc_Plot.curve(time, a_12_sim / char_acc, zorder=1,
                      label=r"$1-2$ Simulation")
    oscAcc_Plot.curve(time, a_32_sim / char_acc, zorder=2,
                      label=r"$3-2$ Simulation")
    oscAcc_Plot.curve(
        time, (a_12_con_mag) / char_acc, thin_curve=True, zorder=3, label=r"$1-2$ Constraint")
    oscAcc_Plot.curve(
        time, (a_32_con_mag) / char_acc, thin_curve=True, zorder=4, label=r"$3-2$ Constraint")
    # Add legend
    oscAcc_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscAcc_Plot.save_plot()

    numLines = 2
    oscAccErr_Plot = PlotStyling(numLines,
                                 r"$t/\tau$", r"Error $(\Delta \dot{\mathbf{U}} \cdot \mathbf{q}) / (U_0 \, \omega)$",
                                 title=None, loglog=False,
                                 outputDir=output_dir, figName="osc-acc-err", eps=epsOutput,
                                 continuousColors=False)

    # Show numerical data points
    oscAccErr_Plot.make_plot()
    oscAccErr_Plot.curve(time, a_12_sim / char_acc -
                         (a_12_con_mag) / char_acc, zorder=1, label=r"$1-2$")
    oscAccErr_Plot.curve(time, a_32_sim / char_acc -
                         (a_32_con_mag) / char_acc, zorder=2, label=r"$3-2$")
    oscAccErr_Plot.set_yaxis_scientific()
    # Add legend
    oscAccErr_Plot.legend(title=r"$Z_0/a =$" + "{}".format(
        fmt(np.max(Z_height))),
        loc='best', ncol=2, bbox_to_anchor=(0.0, 1.0, 0.9, 0.1))
    oscAccErr_Plot.save_plot()

# !SECTION (Plots)


# SECTION: Output data

    # Output data to log file
    file = output_dir + "output.txt"
    with open(file, "w") as f:
        print(f"R_avg = {R_avg}", file=f)
        print(f"Z_height = {Z_height}", file=f)
        print(f"Delta R_2 = {DR_loc[:,-1]}", file=f)
        print("", file=f)

        print("Final frame:", file=f)
        print(f"cos(Delta theta) = {cos_theta[-1]}", file=f)
        print(f"Delta theta = {theta[-1]}", file=f)

        print("Energy:", file=f)
        print(f"E_simple/E0 = {E_simple / E0}", file=f)
        print(f"E_total/E0 = {E_total / E0}", file=f)

# !SECTION (Output data)


# SECTION: For use when being called from command line
if __name__ == "__main__":
    """Main method
    """

    # parse user input
    options, remainder = parser.parse_args(sys.argv[1:])
    relative_path = str(options.u_relative_path)
    output_dir = str(options.u_output_dir)

    # generate plots
    aggregate_plots(relative_path, output_dir)

# !SECTION (For use when being called from command line)
