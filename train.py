# SuperFastPython.com
# example of using starmap() with the process pool
from random import random
from time import sleep
from subprocess import Popen, PIPE
from multiprocessing.pool import Pool
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import mrcfile
import numpy as np
from Docking import gnina_docking, filter_docked_poses_by_correlation
from utils import (
    compute_density_map_in_chimeraX,
    split_sdf_file_to_pdbs,
    group_conformers_to_single_file,
    create_folder,
    delete_extension_from_filename,
    rescale_density_map,
    align_density_maps,
    read_density_data_mrc,
)


# task executed in a worker process
# def task(identifier, value):
#     # report a message
#     print(f'Task {identifier} executing with {value}', flush=True)
#     p = Popen(
#         [
#             f"esdfho{identifier}",
#             str(identifier)
#         ],
#         shell=True,
#         stderr=PIPE,
#     )
#     # block for a moment
#     sleep(value)

#     _, stderr = p.communicate()
#     print(f'Task {identifier} with {value} finished! Return code: {p.returncode}, Stderr: {stderr}', flush=True)
#     # return the generated value
#     return (identifier, value)

# # protect the entry point
# if __name__ == '__main__':
#     # create and configure the process pool
#     with Pool(4) as pool:
#         # prepare arguments
#         items = [(i, random()) for i in range(10)]
#         # execute tasks and process results in order
#         result = pool.starmap(task, items)

#         print(result)
#         for r in result:
#             print(f'Got result: {r}', flush=True)
# process pool is closed automatically


# # print(os.getcwd())

# print(os.path.dirname(os.path.realpath(__file__)))


if __name__ == "__main__":
    # data_dir = os.path.sep + os.path.sep.join(["mnt", "cephfs", "projects", "2023110101_Ligand_fitting_to_EM_maps", "PDBbind", "PDBBind_Zenodo_6408497"])
    # cid = "4rr6"
    # gnina_docking(data_dir, cid)
    # a = np.zeros(3)
    # print(a)
    # with Pool(4) as pool:
    #     # prepare arguments
    #     items = [(i, random()) for i in range(40)]
    #     # execute tasks and process results in order
    #     result = pool.starmap(task, items)

    #     print(result)
    #     for r in result:
    #         print(f'Got result: {r}', flush=True)
    # print("Path:")
    # print(os.path.join(os.getcwd(), "asd", "1.txt"))
    # a = "   ADSFDFdasd    "
    # print(a.strip().lower()
    is_log = True
    chimeraX_log_path = os.path.join(os.getcwd(), "cryoEM_maps", "chimeraX_logs")
    create_folder(chimeraX_log_path)
    script_path = os.path.join(os.getcwd(), "chimeraX_scripts")

    input_ligand = "4ui8_ligand.pdb"
    input_protein = "4ui8_protein_processed.pdb"
    raw_data = os.path.join(os.getcwd(), "cryoEM_maps", "raw_molecule_data")
    ligand_path_full = os.path.join(raw_data, input_ligand)
    protein_path_full = os.path.join(raw_data, input_protein)

    temp_data = os.path.join(
        os.getcwd(),
        "cryoEM_maps",
        f"temp_data_{delete_extension_from_filename(input_ligand)}",
    )
    create_folder(temp_data)

    density_resolution = 3.5
    n_box = 16
    first_density_fname = (
        f"first_dens_{delete_extension_from_filename(input_ligand)}.mrc"
    )
    first_density_path_full = os.path.join(temp_data, first_density_fname)
    p = compute_density_map_in_chimeraX(
        ligand_path_full,
        first_density_path_full,
        density_resolution=density_resolution,
        n_box=n_box,
        is_log=is_log,
        log_path=chimeraX_log_path,
        script_path=script_path,
    )

    _, stderr = p.communicate()

    if p.returncode != 0 or stderr:
        raise RuntimeError(f"Failed to compute first density map for docking: {stderr}")

    conformers_path = os.path.join(temp_data, "raw_conformers_data")
    create_folder(conformers_path)
    sdf_conformers_filename = (
        f"docked_confs_{delete_extension_from_filename(input_ligand)}.sdf"
    )
    sdf_conformers_path_full = os.path.join(conformers_path, sdf_conformers_filename)
    n_pos = 10
    box_extension = 5.0

    path_to_gnina = os.path.join(os.getcwd(), "gnina")
    p = gnina_docking(
        ligand_path_full,
        protein_path_full,
        first_density_path_full,
        sdf_conformers_path_full,
        n_pos,
        box_extension=box_extension,
        path_to_gnina=path_to_gnina,
    )

    _, stderr = p.communicate()

    print(p.returncode)

    if p.returncode != 0:
        raise RuntimeError(f"Failed to perform docking: {stderr}")

    remove_Hs = False
    n_confs, conf_path_full_list = split_sdf_file_to_pdbs(
        sdf_conformers_path_full,
        input_ligand,
        pdb_path=conformers_path,
        remove_Hs=remove_Hs,
    )
    print(f"n_confs: {n_confs}")
    print(f"conf_path_full_list: {conf_path_full_list}")

    threshold_correlation = 0.6
    not_found_corr_value = 0.0
    chimeraX_output_base_filename = "output.tmp"
    chimeraX_output_path = os.path.join(temp_data, "chimeraX_output")
    create_folder(chimeraX_output_path)
    clear_chimeraX_output = False
    write_corrs_to_file = True
    corrs_path_full = os.path.join(temp_data, "docking_correlations.txt")

    corrs, n_appr, appr_dock_path_full_list = filter_docked_poses_by_correlation(
        conf_path_full_list,
        first_density_path_full,
        threshold_correlation=threshold_correlation,
        not_found_corr_value=not_found_corr_value,
        chimeraX_output_base_filename=chimeraX_output_base_filename,
        chimeraX_output_path=chimeraX_output_path,
        density_resolution=density_resolution,
        n_box=n_box,
        is_log=is_log,
        log_path=chimeraX_log_path,
        script_path=script_path,
        clear_chimeraX_output=clear_chimeraX_output,
        write_corrs_to_file=write_corrs_to_file,
        corrs_path_full=corrs_path_full,
    )
    print(f"corrs: {corrs}")
    print(f"n_appr: {n_appr}")
    print(f"appr_dock_path_full_list: {appr_dock_path_full_list}")

    if n_appr == 0:
        raise RuntimeError(
            "No appropariate conformers found based on correlation threshold!"
        )

    conformers_filename = (
        f"group_confs_{delete_extension_from_filename(input_ligand)}.pdb"
    )
    conforemrs_path_full = os.path.join(temp_data, conformers_filename)
    group_conformers_to_single_file(
        appr_dock_path_full_list, conforemrs_path_full, delete_input=False
    )

    final_density_fname = (
        f"final_dens_{delete_extension_from_filename(input_ligand)}.mrc"
    )
    final_density_path_full = os.path.join(temp_data, final_density_fname)
    p = compute_density_map_in_chimeraX(
        conforemrs_path_full,
        final_density_path_full,
        density_resolution=density_resolution,
        n_box=n_box,
        is_log=is_log,
        log_path=chimeraX_log_path,
        script_path=script_path,
    )

    _, stderr = p.communicate()

    if p.returncode != 0 or stderr:
        raise RuntimeError(f"Failed to compute final density map: {stderr}")

    dens, header, voxel_size = read_density_data_mrc(final_density_path_full)
    print(header.origin)
    print(header.nx)
    print(header.ny)
    print(header.nz)
    print(voxel_size)
