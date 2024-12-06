import sys
import getopt
import os
from datetime import datetime, timezone

# append the repo's folder to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    log,
    delete_extension_from_filename,
    extract_filename_from_full_path,
) 
from chimeraX_scripts.chimeraX_utilities import run_chimeraX_command as run_com


"""
This script is used to compute correlation between density map of the given molecule
and target density map inside Chimera. The script accepts the following input arguments:
-i: full path to the input molecule file
-r: resolution for the molecule map generation
-t: full path to the target density map
-l: if provided, contains path to the folder where log file will be written
"""

opts, args = getopt.getopt(sys.argv[1:], "i:r:t:l:")

molecule_path_full = None  # full path to the input molecule file (including its name)
density_resolution = None  # density map resolution for the molecule (in Angrstrom)
target_density_path_full = None # full path to the target density file (including its name)
is_log = False # whether we should write logs for Chimera script
log_path = None # path to the log folder (excluding log file name)
log_fname = None # name of the log file
for opt, arg in opts:
    if opt == "-i":
        molecule_path_full = arg
    elif opt == "-r":
        density_resolution = arg
    elif opt == "-t":
        target_density_path_full = arg
    elif opt == "-l":
        is_log = True
        log_path = arg

# check if the arguments are not None after arguments parsing
assert (
    molecule_path_full
), "Path to the input molecule file is None after argument's parsing."

assert density_resolution, "Density resolution is None after argument's parsing."
assert target_density_path_full, "Path to the input target density file is None after argument's parsing."

if is_log:
    assert log_path, "Path to the log folder is None after argument's parsing."

# create log file name if required
if is_log:
    # extract name of the molecule file
    molecule_fname = extract_filename_from_full_path(molecule_path_full)
    log_fname = (
        datetime.utcnow().strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
        + "_mol_corr_"
        + delete_extension_from_filename(molecule_fname)
        + "_log.txt"
    )

# run Chimera commands
if is_log:
    log("Started Chimera commands.", status="INFO", log_path=log_path, log_filename=log_fname)

run_com(
    "open " + molecule_path_full, is_log=is_log, log_path=log_path, log_filename=log_fname
)  # open molecule file

# NOTE: figure out how to move specify grid sppacing in molmap
run_com(
    "molmap #0 {} gridSpacing 1.0".format(density_resolution),
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  # generate density map for the molecule

run_com(
    "open " + target_density_path_full, is_log=is_log, log_path=log_path, log_filename=log_fname
)  # open target density map

# measure correlation between molecule density map and target density map
# NOTE: the output correlation will be sent to the ReplyLog inside Chimera 
# which we catch by inside the subprocess that runs Chimera's script
run_com(
    "measure correlation #0.1 #1", is_log=is_log, log_path=log_path, log_filename=log_fname
) 

# if all commands run successfully, log info message and stop Chimera
if is_log:
    log(
        "Successfully finished Chimera commands! Exiting...",
        status="INFO",
        log_path=log_path,
        log_filename=log_fname,
    )

run_com(
    "stop now",
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  