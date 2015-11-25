import numpy as np
import qg_transform as qgt

def wave_activity_op(pvg, num_levels, *fields):
    """Wave activity operator. The area of each grid is 1
    
    Args:
        pvg: 2d (potential vorticity) field in physical space. Shape (lat, lon)
        num_levels: integer, number of contour levels
        *fields (optional): variable number of other fields
    
    Returns:
        pvg_mean: Lagrangian-mean pvg
        pvg_A: wave activity for pvg
        *As (optional): wave activitiesfor `*fields1`
    
    See:
        Finite-Amplitude Wave Activity and Diffusive Flux of Potential Vorticity
        in Eddy–Mean Flow Interaction, Nakamura and Zhu, 2010
    """
    pvg1d = pvg.flatten().copy()
    num_points = pvg1d.size
    dn = num_points//num_levels
    n0 = num_points - dn*num_levels
    num_lats, num_lons = pvg.shape
    
    idx = np.argsort(pvg1d)
    pvg1d_sorted = pvg1d[idx]
    del pvg1d
    
    int_pv_contour = np.zeros(num_levels)
    int_pv_lateq   = np.zeros(num_levels)
    int_fields_contour = []
    int_fields_lateq   = []
    fields_sorted = []
    for field in fields:
        fields_sorted.append(field.flatten()[idx])
        int_fields_contour.append(np.zeros(num_levels))
        int_fields_lateq.append(np.zeros(num_levels))

    last_lat_eq = 0
    last_points_within = 0
    pvg_mean = np.zeros(num_levels)
    for i in range(num_levels):
        points_within = n0 + (i+1)*dn
        lat_eq = points_within//num_lons
        pvg_mean[i] = pvg1d_sorted[points_within-1]
        
        int_pv_contour[i] = pvg1d_sorted[last_points_within:points_within].sum()
        int_pv_lateq[i]   = pvg[last_lat_eq:lat_eq].sum()
        for j in range(len(fields)):
            int_fields_contour[j][i] = fields_sorted[j][last_points_within:points_within].sum()
            int_fields_lateq[j][i] = fields[j][last_lat_eq:lat_eq].sum()
        last_lat_eq = lat_eq
        last_points_within = points_within
    pvg_A = -np.cumsum(int_pv_contour - int_pv_lateq)
    As = []
    for i in range(len(fields)):
        As.append(-np.cumsum(int_fields_contour[i] - int_fields_lateq[i]))
    if As:
        return pvg_mean, pvg_A, As
    else:
        return pvg_mean, pvg_A

class TwoLayerWaveActivity(object):
    """analysis finite wave activity budget in a two layer model"""
    def __init__(self, psic, para_dict):
        """
        Args:
            psic: stream function spectrum with shape (time, ky, kx, z)
            para_dict: a dictionary of parameters with keys. For example
                    'beta': 
                    'F':
                    'bot_drag':
                    'U': 1
                    'dt_tune': 0.3
                    'filter_exp': 4
                    'dealiasing':'orszag'
                    'filter_exp': 4
                    'filter_tune':1
        """
        kmax  = psic.shape[1] - 1
        self._psic = psic
        self._ntime = psic.shape[0]
        if psic.shape[3] != 2:
            raise ValueError('Only work for 2-layer model')

        default_para = {'U': 1,
                    'dt_tune'    : 0.3,
                    'filter_exp' : 4,
                    'dealiasing' :'orszag',
                    'filter_tune': 1.0,
                    'filter_type': 'hyperviscous',
                    'Nexp'       : 1.0,
                    'kmax'       : kmax}
        self.params = default_para
        for k in para_dict:
            self.params[k] = para_dict[k]
            
        pvc = qgt.get_PV(psic, self.params['F'])
        self._pvg = qgt.spec2grid(pvc)
        dt = qgt.time_step(pvc, self.params['dt_tune'], self.params['beta'],
                            kmax)
        dt.shape = dt.shape + (1,1,1)
        filter_rate2d = qgt.hypervis_filter_rate(kmax, 1.0, 
            self.params['filter_tune'], self.params['filter_exp'], 
            self.params['dealiasing'], self.params['filter_type'],
            self.params['Nexp'])
        filter_rate2d.shape = (1,) + filter_rate2d.shape + (1,)
        hyp_visc = pvc*filter_rate2d/dt
        self._hyp_visg = qgt.spec2grid(hyp_visc)
        
        self._num_lats = self._pvg.shape[1]
        self._UY = 4*F*np.linspace(-np.pi, np.pi, self._num_lats)[:,np.newaxis]
        self._BETAY = beta*np.linspace(-np.pi, np.pi, self._num_lats)[:,np.newaxis]
        self._eqlats = self._num_lats
        
    def cal_upperlayer(self):
        lang_pv1   = np.zeros(self._eqlats)
        A_pv1      = np.zeros(self._eqlats)
        hypvis_pv1 = np.zeros(self._eqlats)
        for t in range(self._ntime):
            curr_pv1, curr_A, curr_hyvis = wave_activity_op(
                self._pvg[t,...,0]+self._UY+self._BETAY, self._eqlats, 
                self._hyp_visg[t,...,0])
            lang_pv1   += curr_pv1
            A_pv1      += curr_A
            hypvis_pv1 += curr_hyvis[0]
        self.lang_pv1   = lang_pv1/self._ntime
        self.A_pv1      = A_pv1/self._ntime
        self.hypvis_pv1 = hypvis_pv1/self._ntime
        
    def cal_lowerlayer(self):
        lang_pv2   = np.zeros(self._eqlats)
        A_pv2      = np.zeros(self._eqlats)
        hypvis_pv2 = np.zeros(self._eqlats)
        drag_pv2   = np.zeros(self._eqlats)
        
        vorc_lower = qgt.laplacian(self._psic[...,1])
        vorg_drag  = -self.params['bot_drag']*qgt.spec2grid(vorc_lower)
        for t in range(self._ntime):
            curr_pv2, curr_A, budgets = wave_activity_op(
                -(self._pvg[t,...,1]-self._UY+self._BETAY), self._eqlats, 
                self._hyp_visg[t,...,1], vorg_drag[t])
            lang_pv2   += curr_pv2
            A_pv2      += curr_A
            hypvis_pv2 += budgets[0]
            drag_pv2   += budgets[1]
        self.lang_pv2   = -lang_pv2/self._ntime
        self.A_pv2      = A_pv2/self._ntime
        self.hypvis_pv2 = hypvis_pv2/self._ntime
        self.drag_pv2   = drag_pv2/self._ntime
        