"""A collection of functions to do spectral diagnostics in a barotropic flow"""
import qg_transform as qgt


def nonlinear_transfer(psi, S):
    """Evaluate Re(S_kl' * J_kl(psi, S)), ' denote complex conj
    if S is vors, then the result is enstrophy transfer into wavenumber k
    if S is tracers, then the result is tracer variance cascade into k
    
    Inputs:
        psi: spectrum of streamfunction with (t(optional), ky, kx, z(optional))
        S: same shape as psi
    
    Output:
        ks: 1 to kmax (511 in most of my simulations)
        trans: 1d array
    """
    Jg_psi_S = qgt.jacobian(psi, S)
    Js_psi_S = qgt.grid2spec(Jg_psi_S)
    return qgt.prod_spectrum(S, Js_psi_S)

def eddy_mean_transfer(psi, S):
    """Decompose nonlinear transfer term 
           S * J(psi, S)
    into eddy-eddy interaction and eddy-meanflow interaction
    eddy-eddy: S'*J(psi', S')
    eddy-meanflow: S'*(J(psi,S)-J(psi',S')
    in the spectral space. ' means derivation from zonal average
    
    Inputs:
        psi: spectrum of streamfunction with (t(optional), ky, kx, z(optional))
        S: same shape as psi
    
    Output:
        ks: 1 to kmax (511 in most of my simulations)
        eddy_eddy_trans: 1d array
        eddy_mean_trans: 1d array
        total_trans: 1d array, spectral transfer of the whole flow
    """
    psi_eddy = qgt.filter(psi, None, None, True)
    S_eddy   = qgt.filter(S, None, None, True)
    Jg_total = qgt.jacobian(psi, S)
    Jg_eddy  = qgt.jacobian(psi_eddy, S_eddy)
    Jgs_total = qgt.grid2spec(Jg_total)
    Jgs_eddy  = qgt.grid2spec(Jg_eddy)
    k, eddy_eddy_trans = qgt.prod_spectrum(S_eddy, Jgs_eddy)
    k, eddy_mean_trans = qgt.prod_spectrum(S_eddy, Jgs_total-Jgs_eddy)
    k, total_trans = qgt.prod_spectrum(S, Jgs_total)
    return k, eddy_eddy_trans, eddy_mean_trans, total_trans