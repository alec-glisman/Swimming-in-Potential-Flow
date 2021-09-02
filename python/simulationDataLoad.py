# External Dependencies
import os                          # Access system file-tree
import sys                         # Modify system parameters
import tarfile                     # data I/O
from pathlib import Path           # removing filename extension
import warnings

import glob                        # get files matching pattern
import re                          # regex

import numpy as np                 # Data structures
import pandas as pd                # Data structures

# Global options for external dependencies
np.set_printoptions(threshold=sys.maxsize, precision=16)

# Internal Dependencies
from GSDUtil import GSDUtil  # noqa: E402


def find_compressed_data(relative_path_base_dir, sim_parameters_varied):

    data_path = []

    for r in relative_path_base_dir:

        file_paths = glob.glob(r + '/*')  # look for all files in r
        desired_file_paths = []

        for p in sim_parameters_varied:

            # regex matching for tar.xz files of parameter
            regex_str = r + '/.*' + p + '.tar.xz'
            regex = re.compile(fr"{regex_str}")
            matched = list(filter(regex.search, file_paths))

            if (len(matched) == 0):
                warnings.warn(
                    f"Failure to find data. No data path found in parameter: {p}, relative path: {r}")

            desired_file_paths.append(matched)

        # flatten list and output to data_path
        desired_file_paths_flattened = [
            item for sublist in desired_file_paths for item in sublist]
        data_path.append(desired_file_paths_flattened)

    return data_path


def load_compressed_data(tar_data_path):

    # get filename without extension
    filename = Path(tar_data_path)
    filename_wo_ext = filename.with_suffix('').with_suffix('')

    # Decompress data
    try:
        tar_data = tarfile.open(str(tar_data_path))
    except:
        raise IOError(
            f"Input tar_data_path is not a path to tarfile: {tar_data_path}")

    tar_data.extractall(filename_wo_ext)

    # Load data
    gsd_files = []

    try:
        # Loop through all subdirectories in the 'data' directory
        for root, dirs, files in os.walk(filename_wo_ext, topdown=True):
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
            f"Failure to load data. No files found in relPath {filename_wo_ext}")

    # Remove decompressed data
    os.system(f"rm -r {filename_wo_ext}")

    # Return data
    return gsd_files


def parse_loaded_data(gsd_files):

    # NOTE: Assuming each GSD has same number of frames and particles
    nframes = gsd_files[0].trajectory.file.nframes
    nparticles = gsd_files[0].snapshot.particles.N

    # Vector length parameters
    num_input_files = len(gsd_files)
    num_spatial_dimensions = 3
    num_frames = nframes - 2
    num_particles = nparticles

    # Simulation tags
    gsd_paths = []

    # ANCHOR: Load GSD data into np structures
    CoM_disp_comp = np.zeros(
        (num_input_files, num_spatial_dimensions), dtype=np.double)
    CoM_disp = np.zeros(num_input_files, dtype=np.double)

    R_avg = np.zeros_like(CoM_disp, dtype=np.double)
    phaseShift = np.zeros_like(CoM_disp, dtype=np.double)
    U0 = np.zeros_like(CoM_disp, dtype=np.double)
    omega = np.zeros_like(CoM_disp, dtype=np.double)
    epsilon = np.zeros_like(CoM_disp, dtype=np.double)

    # Orientation vector (not unit norm)
    q = np.zeros((num_input_files, num_spatial_dimensions,
                 num_frames), dtype=np.double)
    q0 = np.zeros((num_input_files, num_spatial_dimensions), dtype=np.double)

    # temporal data
    time = np.zeros((num_input_files, num_frames), dtype=np.double)

    # kinematic data (simulation_number, particle_number, dimension_number, frame_number)
    positions = np.zeros((num_input_files, num_particles,
                         num_spatial_dimensions, num_frames), dtype=np.double)
    velocities = np.zeros_like(positions, dtype=np.double)
    accelerations = np.zeros_like(positions, dtype=np.double)

    # orientational displacement
    theta = np.zeros((num_input_files, num_frames), dtype=np.double)
    theta_dot = np.zeros_like(theta, dtype=np.double)
    theta_ddot = np.zeros_like(theta, dtype=np.double)

    for i in range(len(gsd_files)):

        gsd_current = gsd_files[i]
        gsd_paths.append(gsd_current.gsdPath)

        # Data from final frame
        gsd_files[i].snapshot = gsd_current.trajectory.read_frame(nframes - 1)
        CoM_disp_comp[i, :] = gsd_current.snapshot.log['particles/double_position'][1]

        # Data from initial frame (not 0)
        gsd_files[i].snapshot = gsd_current.trajectory.read_frame(1)
        R_avg[i] = gsd_current.snapshot.log['swimmer/R_avg']
        phaseShift[i] = gsd_current.snapshot.log['swimmer/phase_shift']
        U0[i] = gsd_current.snapshot.log['swimmer/U0']
        omega[i] = gsd_current.snapshot.log['swimmer/omega']
        epsilon[i] = U0[i] / R_avg[i] / omega[i]
        CoM_disp_comp[i, :] -= gsd_current.snapshot.log['particles/double_position'][1]

        # Data using both initial and final frame
        CoM_disp[i] = np.linalg.norm(CoM_disp_comp[i, :])

        # Data calculated at each frame
        q_current = np.zeros((3, nframes - 2), dtype=np.double)

        for j in range(1, nframes - 1):  # loop over all frames in GSD file

            snapshot_current = gsd_current.trajectory.read_frame(j)
            jj = j - 1

            q_current[:, jj] = snapshot_current.log['particles/double_position'][0]
            q_current[:, jj] -= snapshot_current.log['particles/double_position'][2]

            time[i, jj] = snapshot_current.log['integrator/t']

            positions[i, :, :,
                      jj] = snapshot_current.log['particles/double_position']
            velocities[i, :, :,
                       jj] = snapshot_current.log['particles/double_velocity']
            accelerations[i, :, :,
                          jj] = snapshot_current.log['particles/double_moment_inertia']

        q[i, :, :] = q_current / np.linalg.norm(q_current, axis=0)
        q0[i, :] = q[i, :, 0]

        # swimmer angular rotation from initial configuration
        cos_theta = np.dot(q[i, :, :].T, q0[i, :])
        theta[i, :] = np.arccos(cos_theta)
        theta_dot[i, :] = np.gradient(theta[i, :], time[i, :])
        theta_ddot[i, :] = np.gradient(theta_dot[i, :], time[i, :])

    # Create dictionary of data
    dict = {
        "gsd_path": gsd_paths,

        "CoM_disp": CoM_disp,
        "CoM_disp_comp": list(CoM_disp_comp),

        "R_avg": R_avg,
        "phase_shift": phaseShift,
        "U0": U0,
        "omega": omega,
        "epsilon": epsilon,

        "time": list(time),
        "positions": list(positions),
        "velocities": list(velocities),
        "accelerations": list(accelerations),

        "initial_time": list(time[:, 0]),
        "initial_positions": list(positions[:, :, :, 0]),
        "initial_velocities": list(velocities[:, :, :, 0]),
        "initial_accelerations": list(accelerations[:, :, :, 0]),

        "final_time": list(time[:, -1]),
        "final_positions": list(positions[:, :, :, -1]),
        "final_velocities": list(velocities[:, :, :, -1]),
        "final_accelerations": list(accelerations[:, :, :, -1]),

        "q": list(q),
        "q0": list(q0),

        "theta": list(theta),
        "theta_dot": list(theta_dot),
        "theta_ddot": list(theta_ddot),

        "final_theta": list(theta[:, -1])
    }

    # Output data
    df = pd.DataFrame(data=dict)
    return df


def gsd_df(relative_path_base_dir, sim_parameters_varied):

    data_df = []

    # Find path to relevant compressed data
    comp_data_path = find_compressed_data(
        relative_path_base_dir, sim_parameters_varied)

    for i in range(len(comp_data_path)):  # loop over data source directories

        data_source_gsd_files = []

        for j in range(len(comp_data_path[i])):  # loop over parameters varied

            # Load GSD files from compressed data
            gsd_current = load_compressed_data(comp_data_path[i][j])

            # Parse GSD files
            if (len(gsd_current) > 0):
                data_source_gsd_files.append(gsd_current)

        # Flatten output
        data_source_gsd_files = [
            item for sublist in data_source_gsd_files for item in sublist]

        # Parse GSD files
        data_df.append(parse_loaded_data(data_source_gsd_files))

    # Combine all data
    all_data_df = pd.concat(data_df)
    return all_data_df
