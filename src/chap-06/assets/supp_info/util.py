import os
import re
import numpy as np
import scipy.constants
import pandas as pd


BOHR = scipy.constants.physical_constants["Bohr radius"][0] / 1e-10


def read_mat(lines):
    if isinstance(lines, str):
        lines = lines.split("\n")
    return np.array([re.sub("[\[\]]", " ", l).split() for l in lines], dtype=float)


def read_comp(mat, atol=1e-5, rtol=1e-4):
    assert np.allclose(mat, mat.T, atol=atol, rtol=rtol)
    assert mat.shape == (3, 3)
    return (0.5 * mat + 0.5 * mat.T)[(0, 1, 2, 0, 1, 2), (0, 1, 2, 1, 2, 0)]


def get_iso(arr):
    if arr.shape == (3, 3):
        arr = read_comp(arr)
    assert arr.shape == (6, )
    return 1 / 3 * arr[:3].sum()


def get_aniso(arr, allow_diag=False):
    if arr.shape == (3, 3):
        arr = read_comp(arr)
    assert arr.shape in [(6, ), (3, )]
    if arr.shape == (6, ):
        xx, yy, zz, xy, yz, zx = arr
    else:
        assert allow_diag, "Anisotropic polarizability should include xy, yz, zx components."
        xx, yy, zz = arr
        xy, yz, zx = 0, 0, 0
    return np.sqrt(0.5) * np.sqrt((xx - yy)**2 + (yy - zz)**2 + (zz - xx)**2 + 6 * (xy**2 + yz**2 + zx**2))


def read_by_prompt(lines, token):
    for n, l in enumerate(lines):
        if l.startswith(token):
            return read_mat(lines[n+1:n+4])
    assert False


def get_df_err(df, df_sub, ref=None):
    if ref is None:
        ref = df_sub
    # In case if column headers of reference is not the same to df
    if isinstance(df_sub, pd.DataFrame):
        df_sub.columns = df.columns
        ref.columns = df.columns
    df_err_val = df.sub(df_sub, axis="index")
    df_err_rel = df_err_val.div(ref, axis="index") * 100
    df_z = df - df.mean()
    ref_z = df_sub - df_sub.mean()
    df_err = {
        "MaxE/A^3": df_err_val.abs().max(),
        "MAD/A^3": df_err_val.abs().mean(),
        "RMSD/A^3": (df_err_val**2).mean()**0.5,
        "RelMaxE/%": df_err_rel.abs().max(),
        "RelMAD/%": df_err_rel.abs().mean(),
        "RelRMSD/%": (df_err_rel**2).mean()**0.5,
    }
    # In case df is pd.Series instead of pd.DataFrame 
    try:
        return pd.DataFrame(df_err).T
    except ValueError:
        return pd.Series(df_err)

    
def get_rmsre_3comp(df):
    assert list(df.columns) == ["xx", "yy", "zz"]
    return np.sqrt((df.loc["RelRMSD/%"].mean()**2))


def get_relrmsd_3comp(df):
    assert list(df.columns) == ["xx", "yy", "zz"]
    return np.sqrt((df.loc["RelRMSD/%"]**2).mean())
