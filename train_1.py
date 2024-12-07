import os

import scipy.ndimage
from utils import (
    compute_mol_map_correlation_in_chimera,
    delete_extension_from_filename,
    read_density_data_mrc,
    calculate_center_of_mass_of_density,
)
from subprocess import Popen, PIPE
import mrcfile
import numpy as np
import scipy

mol_path_full = os.path.join(
    os.getcwd(), "cryoEM_maps", "raw_molecule_data", "4ui8_ligand.pdb"
)
target_path_full = os.path.join(
    os.getcwd(), "cryoEM_maps", "temp_data_4ui8_ligand", "first_dens_4ui8_ligand.mrc"
)
log_path = os.path.join(os.getcwd(), "cryoEM_maps", "chimera_logs")
output_path_full = os.path.join(
    os.getcwd(), "cryoEM_maps", "temp_data_4ui8_ligand", "chimera_output.txt"
)

# p = compute_mol_map_correlation_in_chimera(
#     mol_path_full,
#     target_path_full,
#     is_log=True,
#     log_path=log_path,
#     output_path_full=output_path_full
# )

# _, stderr = p.communicate()
# print(p.returncode)
# print(stderr)

# script_line = os.path.join(os.getcwd(), "chimera_scripts", "chimera_mol_density_corr.py")
# command = [
#         "chimera",
#         "--nogui",
#         "--nostatus",
#         "--script",
#         f"{script_line} -i {mol_path_full} -r {1.0} -t {target_path_full} -l {log_path}"
#         ]

# file = open(output_path_full, "w")
# p = Popen(
#     command,
#     stdout=file,
#     stderr=PIPE
# )


# stdout, stderr = p.communicate()

# file.close()

# print(p.returncode)
# print(stdout)
# print(stderr)

# input_ligand = "4ui8_ligand.pdb"
# temp_data = os.path.join(os.getcwd(), "cryoEM_maps", f"temp_data_{delete_extension_from_filename(input_ligand)}")
# scaled_density_fname = f"scaled_dens_{delete_extension_from_filename(input_ligand)}.mrc"
# density_fname = f"dens_{delete_extension_from_filename(input_ligand)}.mrc"
# scaled_density_path_full = os.path.join(temp_data, scaled_density_fname)
# density_path_full = os.path.join(temp_data, density_fname)

# # with mrcfile.open(scaled_density_path_full, 'r+') as mrc1:
# #     with mrcfile.open(density_path_full, 'r+') as mrc:
# #         mrc1.header.origin = mrc.header.origin
# #         print(mrc.header.origin)
# #         print(mrc1.header.origin)


# map_1 = os.path.join(os.getcwd(), "cryoEM_maps", "temp_data_4ui8_ligand", "dens_4ui8_ligand.mrc")
# map_2 = os.path.join(os.getcwd(), "cryoEM_maps", "temp_data_4ui8_ligand", "scaled_dens_4ui8_ligand.mrc")
# dnesity_1, header_1, voxel_1 = read_density_data_mrc(map_1)
# origin_1 = header_1.origin.tolist()
# voxel_size_1 = voxel_1.x
# print(voxel_1.x)
# print(voxel_1.y)
# print(voxel_1.z)
# print(origin_1)
# print(header_1.nx)
# print(header_1.ny)
# print(header_1.nz)
# center_1 = calculate_center_of_mass_of_density(dnesity_1)
# print(center_1)

# dnesity_2, header_2, voxel_2 = read_density_data_mrc(map_2)
# origin_2 = header_2.origin.tolist()
# voxel_size_2 = voxel_2.x
# print(voxel_2.x)
# print(voxel_2.y)
# print(voxel_2.z)
# print(origin_2)
# print(header_2.nx)
# print(header_2.ny)
# print(header_2.nz)
# center_2 = calculate_center_of_mass_of_density(dnesity_2)
# print(center_2)

# print(dnesity_2.shape)


# non_zero = np.nonzero(dnesity_2)

# print(dnesity_2[non_zero[0][0], non_zero[1][0], non_zero[2][0]])
# print(non_zero[0][0] * voxel_2.x , non_zero[1][0] * voxel_2.y, non_zero[2][0] * voxel_2.z)

# non_zero = np.nonzero(dnesity_1)

# print(dnesity_1[non_zero[0][0], non_zero[1][0], non_zero[2][0]])
# print(non_zero[0][0] * voxel_1.x , non_zero[1][0] * voxel_1.y, non_zero[2][0] * voxel_1.z)

# # print(np.nonzero(a))
# ccc = (dnesity_2 * np.mgrid[0:dnesity_2.shape[0], 0:dnesity_2.shape[1], 0:dnesity_2.shape[2]]).sum(1).sum(1).sum(1) / dnesity_2.sum()

# print(f"New center: {ccc}")


# print(f"New new center: {scipy.ndimage.center_of_mass(dnesity_2)}")
temp_data = os.path.join(os.getcwd(), "cryoEM_maps", "temp_data_4ui8_ligand")
output_dens_full = os.path.join(temp_data, "map.mrc")


dens, header, voxel_size = read_density_data_mrc(output_dens_full)

print(header.origin)
print(header.nx)
print(header.ny)
print(header.nz)
print(voxel_size)
