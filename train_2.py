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
    read_molecule,
    save_all_conformers_to_pdb,
    read_density_data_mrc,
)


def generate_conformers(
    mol,
    n_confs,
    add_Hs=False,
    use_small_ring_torsions=False,
    prune_rms_tresh=1.0,
    random_seed=0xF00D,
):
    """
    Generates n_confs conformers for the give molecule using RDKit library.

    Params:
    mol - input RDKit molecule object
    n_confs - number of conformers to generate
    add_Hs - whether to add hydrogen atoms to the given molecule before conformers generation
    use_small_ring_torsions - whether to include additional small ring torsion potentials
    prune_rms_tresh - threshold value for RMSD pruning of the conformers
    random_seed - seed for conformers generation algorithm (for reproducibility)

    Returns:
    mol - RDKit molecule object with generated conformers
    """

    if add_Hs:  # add hydrogen atoms if needed
        mol = Chem.AddHs(mol)

    # specify parameters for the conformers generation
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = use_small_ring_torsions
    params.pruneRMsThresh = prune_rms_tresh
    params.randomSeed = random_seed
    params.useRandomCoords = False

    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    return mol


if __name__ == "__main__":
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
        f"RDKit_temp_data_{delete_extension_from_filename(input_ligand)}",
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

    mol = read_molecule(ligand_path_full, remove_Hs=False)

    original_coords = mol.GetConformer().GetPositions()
    original_centroid = np.mean(original_coords, axis=0)

    original_mol = Chem.Mol(mol)  # copy original molecule for future allignment

    # generate conformers
    n_confs = 10
    mol = generate_conformers(
        mol,
        n_confs=n_confs,
        add_Hs=False,
        use_small_ring_torsions=False,
        prune_rms_tresh=1.0,
    )

    conformer_ids = [x.GetId() for x in mol.GetConformers()]
    print(conformer_ids)

    for conf_id in conformer_ids:
        conformer = mol.GetConformer(conf_id)
        conformer_coords = conformer.GetPositions()
        conformer_centroid = np.mean(conformer_coords, axis=0)
        shift = original_centroid - conformer_centroid
        new_coords = conformer_coords + shift
        for i, coord in enumerate(new_coords):
            conformer.SetAtomPosition(i, coord)

    # for conf_id in conformer_ids:
    #     conformer = mol.GetConformer(conf_id)
    #     AllChem.AlignMol(mol, original_mol, prbCid=conf_id, refCid=0)  # Align conformer to original molecule

    conformers_path = os.path.join(temp_data, "raw_conformers_data")
    create_folder(conformers_path)
    base_conformer_filename = delete_extension_from_filename(input_ligand) + ".pdb"
    _, conformer_path_list = save_all_conformers_to_pdb(
        mol, base_conformer_filename, pdb_path=conformers_path
    )

    print(conformer_path_list)

    threshold_correlation = 0.6
    not_found_corr_value = 0.0
    chimeraX_output_base_filename = "output.tmp"
    chimeraX_output_path = os.path.join(temp_data, "chimeraX_output")
    create_folder(chimeraX_output_path)
    clear_chimeraX_output = False
    write_corrs_to_file = True
    corrs_path_full = os.path.join(temp_data, "docking_correlations.txt")

    corrs, n_appr, appr_dock_path_full_list = filter_docked_poses_by_correlation(
        conformer_path_list,
        first_density_path_full,
        threshold_correlation=threshold_correlation,
        not_found_corr_value=not_found_corr_value,
        chimeraX_output_base_filename=chimeraX_output_base_filename,
        chimeraX_output_path=chimeraX_output_path,
        density_resolution=density_resolution,
        n_box=n_box,
        is_log=is_log,
        log_path=chimeraX_log_path,
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
