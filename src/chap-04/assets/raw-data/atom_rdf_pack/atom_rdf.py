import numpy as np
from pyscf import gto, data


ngrids = 50000
rad = np.linspace(0, 10 / data.nist.BOHR, ngrids)
coords = np.zeros((ngrids, 3))
coords[:, 0] = rad
factors = 4 * np.pi * rad**2


def get_ao(mol):
    return mol.eval_gto("GTOval_sph_deriv2", coords=coords)


def get_rdf(mol, rdm1):
    ao = get_ao(mol)
    ao_0 = ao[0]
    ao_1 = ao[1:4]
    ao_2 = ao[[4, 7, 9]]
    ao_0_dm = ao_0 @ rdm1
    ao_1_dm = ao_1 @ rdm1
    rho = np.einsum("g, gu, gu -> g", factors, ao_0_dm, ao_0, optimize=True)
    grd = np.einsum("g, gu, tgu -> tg", factors, ao_0_dm, ao_1, optimize=True)
    grd = 2 * np.linalg.norm(grd, axis=0)
    lr_1 = np.einsum("g, tgu, tgu -> g", factors, ao_1_dm, ao_1, optimize=True)
    lr_2 = np.einsum("g, gu, tgu -> g", factors, ao_0_dm, ao_2, optimize=True)
    lr = 2 * (lr_1 + lr_2)
    return {
        "RHO": rho,
        "GRD": grd,
        "LR": lr
    }


def get_rmsd_error(a, b, scale=1):
    if isinstance(scale, str):
        if scale.upper() == "RHO":
            scale = 0.009943368
        elif scale.upper() == "GRD":
            scale = 0.092398036
        elif scale.upper() == "LR":
            scale = 1.445110833
        else:
            assert False
    return 1 / scale * np.sqrt(((a - b)**2).sum() / a.size)
