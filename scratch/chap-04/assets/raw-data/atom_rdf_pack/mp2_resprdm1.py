import numpy
import numpy as np
from functools import reduce
from pyscf import mp
from pyscf import lib
from pyscf.scf import cphf
from pyscf.ao2mo import _ao2mo
from pyscf.grad.mp2 import _shell_prange, _index_frozen_active
from pyscf.mp.mp2 import _gamma1_intermediates


def response_dm1(mymp):
    # general information
    mol = mymp.mol
    t2 = mymp.t2
    with_frozen = not ((mymp.frozen is None)
                       or (isinstance(mymp.frozen, (int, numpy.integer)) and mymp.frozen == 0)
                       or (len(mymp.frozen) == 0))
    OA, VA, OF, VF = _index_frozen_active(mymp.get_frozen_mask(), mymp.mo_occ)
    orbo = mymp.mo_coeff[:,OA]
    orbv = mymp.mo_coeff[:,VA]
    nao, nocc = orbo.shape
    nvir = orbv.shape[1]
    atmlst = range(mol.natm)
    offsetdic = mol.offset_nr_by_atom()
    diagidx = numpy.arange(nao)
    diagidx = diagidx*(diagidx+1)//2 + diagidx
    
    # doo, dvv
    mp.nmo = nocc + nvir
    mp.nocc = nocc
    d1 = _gamma1_intermediates(mp, t2)
    doo, dvv = d1
    
    # part_dm2
    part_dm2 = _ao2mo.nr_e2(t2.reshape(nocc**2,nvir**2),
                            numpy.asarray(orbv.T, order='F'), (0,nao,0,nao),
                            's1', 's1').reshape(nocc,nocc,nao,nao)
    part_dm2 = (part_dm2.transpose(0,2,3,1) * 4 -
                part_dm2.transpose(0,3,2,1) * 2)

    # Imat
    Imat = numpy.zeros((nao,nao))
    max_memory = max(0, mymp.max_memory - lib.current_memory()[0])
    blksize = max(1, int(max_memory*.9e6/8/(nao**3*2.5)))
    for k, ia in enumerate(atmlst):
        shl0, shl1, p0, p1 = offsetdic[ia]
        ip1 = p0
        for b0, b1, nf in _shell_prange(mol, shl0, shl1, blksize):
            ip0, ip1 = ip1, ip1 + nf
            dm2buf = lib.einsum('pi,iqrj->pqrj', orbo[ip0:ip1], part_dm2)
            dm2buf+= lib.einsum('qi,iprj->pqrj', orbo, part_dm2[:,ip0:ip1])
            dm2buf = lib.einsum('pqrj,sj->pqrs', dm2buf, orbo)
            dm2buf = dm2buf + dm2buf.transpose(0,1,3,2)
            dm2buf = lib.pack_tril(dm2buf.reshape(-1,nao,nao)).reshape(nf,nao,-1)
            dm2buf[:,:,diagidx] *= .5

            shls_slice = (b0,b1,0,mol.nbas,0,mol.nbas,0,mol.nbas)
            eri0 = mol.intor('int2e', aosym='s2kl', shls_slice=shls_slice)
            Imat += lib.einsum('ipx,iqx->pq', eri0.reshape(nf,nao,-1), dm2buf)
            eri0 = None
            dm2buf = None


    mo_coeff = mymp.mo_coeff
    mo_energy = mymp._scf.mo_energy
    nao, nmo = mo_coeff.shape
    nocc = numpy.count_nonzero(mymp.mo_occ > 0)
    Imat = reduce(numpy.dot, (mo_coeff.T, Imat, mymp._scf.get_ovlp(), mo_coeff)) * -1

    dm1mo = numpy.zeros((nmo,nmo))
    if with_frozen:
        dco = Imat[OF[:,None],OA] / (mo_energy[OF,None] - mo_energy[OA])
        dfv = Imat[VF[:,None],VA] / (mo_energy[VF,None] - mo_energy[VA])
        dm1mo[OA[:,None],OA] = doo + doo.T
        dm1mo[OF[:,None],OA] = dco
        dm1mo[OA[:,None],OF] = dco.T
        dm1mo[VA[:,None],VA] = dvv + dvv.T
        dm1mo[VF[:,None],VA] = dfv
        dm1mo[VA[:,None],VF] = dfv.T
    else:
        dm1mo[:nocc,:nocc] = doo + doo.T
        dm1mo[nocc:,nocc:] = dvv + dvv.T

    dm1 = reduce(numpy.dot, (mo_coeff, dm1mo, mo_coeff.T))
    vhf = mymp._scf.get_veff(mymp.mol, dm1) * 2
    Xvo = reduce(numpy.dot, (mo_coeff[:,nocc:].T, vhf, mo_coeff[:,:nocc]))
    Xvo+= Imat[:nocc,nocc:].T - Imat[nocc:,:nocc]
    dm1mo += _response_dm1(mymp, Xvo)
    
    dm1 = reduce(numpy.dot, (mo_coeff, dm1mo, mo_coeff.T))
    hf_dm1 = mymp._scf.make_rdm1(mymp.mo_coeff, mymp.mo_occ)
    dm1 += hf_dm1
    dm1 = 0.5 * (dm1 + dm1.T)
    return dm1


def _response_dm1(mymp, Xvo):
    nvir, nocc = Xvo.shape
    nmo = nocc + nvir
    mo_energy = mymp._scf.mo_energy
    mo_occ = mymp.mo_occ
    mo_coeff = mymp.mo_coeff
    def fvind(x):
        x = x.reshape(Xvo.shape)
        dm = reduce(numpy.dot, (mo_coeff[:,nocc:], x, mo_coeff[:,:nocc].T))
        v = mymp._scf.get_veff(mymp.mol, dm + dm.T)
        v = reduce(numpy.dot, (mo_coeff[:,nocc:].T, v, mo_coeff[:,:nocc]))
        return v * 2
    dvo = cphf.solve(fvind, mo_energy, mo_occ, Xvo, max_cycle=30)[0]
    dm1 = numpy.zeros((nmo,nmo))
    dm1[nocc:,:nocc] = dvo
    dm1[:nocc,nocc:] = dvo.T
    return dm1


def get_mp2_resp_rdm1_and_etot(mf_scf):
    mf_mp = mp.MP2(mf_scf).run(with_t2=True)
    rdm1 = response_dm1(mf_mp)
    return rdm1, mf_mp.e_tot
