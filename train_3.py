import os
from itertools import repeat
from utils import (
    compute_density_map_in_chimeraX,
    create_folder,
    compute_mol_map_correlation_in_chimeraX,
    extract_correlation_from_chimera_file,
)

dirname = "/1/2/3/4"
print(os.path.dirname(os.path.dirname(dirname)))

for i in repeat({"a": 1, "b": 2}, 3):
    print(i)
# if return_code_chim != 0 or (err_chim and "Cannot find consistent set of bond" not in err_chim):
#     raise RuntimeError("Chimera's subprocess finished with errors: " + err_chim)
