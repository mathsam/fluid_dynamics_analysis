import numpy as np
import qg_transform as qgt
import scipy.ndimage.interpolation as interp

def wrap_field_NS(field, border_y):
    """Expland a doubly-periodic field in Northern and Sothern borders
    
    Args:
        field: 2d numpy array with shape (num_lats, num_lons)
        border_y: width of the fields that add to NS borders
    
    Returns:
        expand_field: 2d numpy array with shape 
                      (num_lats+2*border_y, num_lons)
    """
    num_lats, num_lons = field.shape
    expand_field = np.zeros((num_lats+2*border_y, num_lons))
    expand_field[:border_y,:] = field[-border_y:,:]
    expand_field[border_y:-border_y,:] = field
    expand_field[-border_y:,:] = field[:border_y,:]
    return expand_field

def wave_activity_op(vor_org, num_levels, fields=None, zoom=4, border_y=128, beta=1.0):
    """Wave activity operator. The domain is [0, 2*pi]x[0, 2*pi]
    
    Args:
        vor_org: 2d (relative vorticity) field in physical space. Shape (lat, lon)
                 potential vorticity is vor_org + beta*y
        num_levels: integer, number of contour levels
        *fields (optional): variable number of other fields
        zoom: integer, the factor of interpolation field to higher resolution
        border_y: add fields of width `border_y` to the Northern and Southern
                  boundaries, and throw away them after done calculation
        beta: float
    
    Returns:
        pvg_mean: Lagrangian-mean pvg
        pvg_A: wave activity for pvg
        fields_mean(optional): Lagrangian-mean along pvg for `*fields`
        As (optional): wave activitiesfor `*fields1`
    
    See:
        Finite-Amplitude Wave Activity and Diffusive Flux of Potential Vorticity
        in Eddy-Mean Flow Interaction, Nakamura and Zhu, 2010
    """
    pvg = interp.zoom(vor_org, zoom, mode='wrap')
    unwrap_num_lats = pvg.shape[0]
    pvg = wrap_field_NS(pvg, border_y)
    num_lats, num_lons = pvg.shape
    num_wrap_levels = int(np.floor(float(border_y)/float(unwrap_num_lats)*num_levels))
    num_levels += num_wrap_levels*2
    dy = 2*np.pi/unwrap_num_lats
    dx = 2*np.pi/pvg.shape[1]
    ys = np.linspace(-border_y*dy, (border_y+unwrap_num_lats-1)*dy,
                     num_lats)
    pvg += (beta*ys)[:,np.newaxis]
    pvg1d = pvg.flatten().copy()
    num_points = pvg1d.size
    dn = num_points//num_levels
    n0 = num_points - dn*num_levels
    
    idx = np.argsort(pvg1d)
    pvg1d_sorted = pvg1d[idx]
    del pvg1d
    
    int_pv_contour = np.zeros(num_levels)
    int_pv_lateq   = np.zeros(num_levels)
    int_fields_contour = []
    int_fields_lateq   = []
    fields_sorted = []
    fields_mean = []
    for i_field in fields:
        field = wrap_field_NS(interp.zoom(i_field, zoom, mode='wrap'), border_y)
        fields_sorted.append(field.flatten()[idx])
        int_fields_contour.append(np.zeros(num_levels))
        int_fields_lateq.append(np.zeros(num_levels))
        fields_mean.append(np.zeros(num_levels))

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
            fields_mean[j][i] = int_fields_contour[j][i]/(points_within-last_points_within)
            int_fields_lateq[j][i] = fields[j][last_lat_eq:lat_eq].sum()
            
        last_lat_eq = lat_eq
        last_points_within = points_within
    pvg_A = -np.cumsum(int_pv_contour - int_pv_lateq)
    As = []
    pvg_mean = pvg_mean[num_wrap_levels:-num_wrap_levels].copy()
    
    frac_dS_Lx = dx*dy/2/np.pi
    pvg_A = pvg_A[num_wrap_levels:-num_wrap_levels].copy()*frac_dS_Lx
    for i in range(len(fields)):
        As.append(
            (-np.cumsum(int_fields_contour[i] - int_fields_lateq[i]))
            [num_wrap_levels:-num_wrap_levels]*frac_dS_Lx)
        fields_mean[i] = fields_mean[i][num_wrap_levels:-num_wrap_levels].copy()
    if As:
        return pvg_mean, pvg_A, fields_mean, As
    else:
        return pvg_mean, pvg_A

def time_mean_wave_activity_op(pvg, num_levels, fields=None, zoom=4, border_y=128, beta=1.0):
    """Wave activity operator. The area of each grid is 1
    First apply wave_activity_op to a snapshot, then time average the snapshots
    
    Args:
        pvg: 3d (potential vorticity) field in physical space. Shape (time, lat, lon)
        num_levels: integer, number of contour levels
        *fields (optional): variable number of other fields
        zoom: integer, the factor of interpolation field to higher resolution
        border_y: add fields of width `border_y` to the Northern and Southern
                  boundaries, and throw away them after done calculation
    
    Returns:
        pvg_mean: Lagrangian-mean pvg
        pvg_A: wave activity for pvg
        fields_mean(optional): Lagrangian-mean along pvg for `*fields`
        As (optional): wave activitiesfor `*fields1`
    
    See:
        Finite-Amplitude Wave Activity and Diffusive Flux of Potential Vorticity
        in Eddy-Mean Flow Interaction, Nakamura and Zhu, 2010
    """
    if pvg.ndim != 3:
        raise TypeError('only work for field with shape (time, lat, lon)')
    num_t = pvg.shape[0]
    if fields is not None:
        num_fields = len(fields)
    else:
        num_fields = 0
    pvg_mean = np.zeros(num_levels)
    pvg_A    = np.zeros(num_levels)
    fields_mean = []
    As          = []
    for i in range(num_fields):
        fields_mean.append(np.zeros(num_levels))
        As.append(np.zeros(num_levels))
    for t in range(num_t):
        curr_fields = []
        for i in range(num_fields):
            curr_fields.append(fields[i][t])
        vals = wave_activity_op(pvg[t], num_levels, curr_fields, 
            zoom, border_y, beta)
        pvg_mean += vals[0]
        pvg_A    += vals[1]
        for i in range(num_fields):
            fields_mean[i] += vals[2][i]
            As[i]          += vals[3][i]
    pvg_mean /= num_t
    pvg_A    /= num_t
    for i in range(num_fields):
        fields_mean[i] /= num_t
        As[i]          /= num_t
    if fields:
        return pvg_mean, pvg_A, fields_mean, As
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
        vorc= qgt.get_vorticity(psic)
        tauc= 0.5*(psic[...,0]-psic[...,1])
        self._pvg = qgt.spec2grid(pvc)
        self._dt = qgt.time_step(pvc, self.params['dt_tune'], 
            self.params['beta'], kmax)
        self._dt.shape = self._dt.shape + (1,1,1)
        filter_rate2d = qgt.hypervis_filter_rate(kmax, 1.0, 
            self.params['filter_tune'], self.params['filter_exp'], 
            self.params['dealiasing'], self.params['filter_type'],
            self.params['Nexp'])
        filter_rate2d.shape = (1,) + filter_rate2d.shape + (1,)
        hyp_visc = pvc*filter_rate2d/self._dt
        self._hyp_visg = qgt.spec2grid(hyp_visc)
        self._hyp_vorg = qgt.spec2grid(vorc*filter_rate2d/self._dt)
        self._hyp_btvorg = qgt.spec2grid(
            np.mean(vorc, -1)*filter_rate2d[...,0]/self._dt[...,0])
        self._hyp_taug = qgt.spec2grid(
            tauc*filter_rate2d[...,0]/self._dt[...,0])
        
        self._num_lats = self._pvg.shape[1]
        self._UY = 4*self.params['F']*np.linspace(-np.pi, np.pi, self._num_lats)[:,np.newaxis]
        self._BETAY = self.params['beta']*np.linspace(-np.pi, np.pi, self._num_lats)[:,np.newaxis]
        self._eqlats = self._num_lats
    
    def cal_barotropic(self):
        bt_vor = np.mean(self._pvg, -1)
        lower_vorg = qgt.spec2grid(qgt.get_vorticity(self._psic[...,1]))
        tauc = 0.5*(self._psic[...,0]-self._psic[...,1])
        minus_J_tau_vortau = -qgt.jacobian(tauc, qgt.laplacian(tauc))
        minus_U_vor_dtaudx = -self.params['U']*qgt.spec2grid(
            qgt.laplacian(qgt.partial_x(tauc)))
        mean_btpv, btpv_A, fs, As = time_mean_wave_activity_op(bt_vor, 
            self._eqlats, 
            [self._hyp_btvorg, minus_J_tau_vortau, 
            -self.params['bot_drag']*lower_vorg/2,
            minus_U_vor_dtaudx], beta=self.params['beta'])
        self.lang_btpv = mean_btpv
        self.A_btpv    = btpv_A
        self.hypvis_btpv    = As[0]
        self.lang_hypvis_btpv   = fs[0]
        self.mJ_tau_vortau  = As[1]
        self.lang_mJ_tau_vortau = fs[1]
        self.drag_btpv      = As[2]
        self.lang_drag_btpv     = fs[2]
        self.mU_vor_dtaudx  = As[3]
        self.lang_mU_vor_dtaudx = fs[3]
        
    def cal_upperlayer(self):
        lang_pv1   = np.zeros(self._eqlats)
        A_pv1      = np.zeros(self._eqlats)
        hypvis_pv1 = np.zeros(self._eqlats)
        hypvis_vor1= np.zeros(self._eqlats)
        hypvis_tau1= np.zeros(self._eqlats)
        for t in range(self._ntime):
            curr_pv1, curr_A, _, curr_hyvis = wave_activity_op(
                self._pvg[t,...,0]+self._UY+self._BETAY, self._eqlats, 
                self._hyp_visg[t,...,0], self._hyp_vorg[t,...,0],
                self._hyp_taug[t]*2*self.params['F'])
            lang_pv1   += curr_pv1
            A_pv1      += curr_A
            hypvis_pv1 += curr_hyvis[0]
            hypvis_vor1+= curr_hyvis[1]
            hypvis_tau1+= curr_hyvis[2]
            
        self.lang_pv1   = lang_pv1/self._ntime
        self.A_pv1      = A_pv1/self._ntime
        self.hypvis_pv1 = hypvis_pv1/self._ntime
        self.hypvis_vor1= hypvis_vor1/self._ntime
        self.hypvis_tau1 = hypvis_tau1/self._ntime
        
    def cal_lowerlayer(self):
        lang_pv2   = np.zeros(self._eqlats)
        A_pv2      = np.zeros(self._eqlats)
        hypvis_pv2 = np.zeros(self._eqlats)
        drag_pv2   = np.zeros(self._eqlats)
        hypvis_vor2= np.zeros(self._eqlats)
        hypvis_tau2= np.zeros(self._eqlats)
        
        vorc_lower = qgt.laplacian(self._psic[...,1])
        vorg_drag  = -self.params['bot_drag']*qgt.spec2grid(vorc_lower)
        for t in range(self._ntime):
            curr_pv2, curr_A, _, budgets = wave_activity_op(
                -(self._pvg[t,...,1]-self._UY+self._BETAY), self._eqlats, 
                self._hyp_visg[t,...,1], vorg_drag[t],
                self._hyp_vorg[t,...,1],
                -self._hyp_taug[t]*2*self.params['F'])
            lang_pv2   += curr_pv2
            A_pv2      += curr_A
            hypvis_pv2 += budgets[0]
            drag_pv2   += budgets[1]
            hypvis_vor2 += budgets[2]
            hypvis_tau2 += budgets[3]
        self.lang_pv2   = -lang_pv2/self._ntime
        self.A_pv2      = A_pv2/self._ntime
        self.hypvis_pv2 = hypvis_pv2/self._ntime
        self.drag_pv2   = drag_pv2/self._ntime
        self.hypvis_vor2 = hypvis_vor2/self._ntime
        self.hypvis_tau2 = hypvis_tau2/self._ntime

class OneLayerWaveActivity(object):
    """analysis finite wave activity budget in a 1-layer barotropic model (F=0)"""
    def __init__(self, psic, para_dict):
        """
        Args:
            psic: stream function spectrum with shape (time, ky, kx)
            para_dict: a dictionary of parameters with keys. For example
                    'beta': 
                    'bot_drag':
                    'dt_tune': 0.3
                    'filter_exp': 4
                    'dealiasing':'orszag'
                    'filter_exp': 4
                    'filter_tune':1
        """
        kmax  = psic.shape[1] - 1
        self._psic = psic
        self._ntime = psic.shape[0]
        if psic.ndim != 3:
            raise ValueError('Only work for 1-layer model')

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
            
        self._vorc= qgt.get_vorticity(psic)
        self._dt = qgt.time_step(self._vorc, self.params['dt_tune'], 
            self.params['beta'], kmax)
        self._dt.shape = self._dt.shape + (1,1)
        filter_rate2d = qgt.hypervis_filter_rate(kmax, 1.0, 
            self.params['filter_tune'], self.params['filter_exp'], 
            self.params['dealiasing'], self.params['filter_type'],
            self.params['Nexp'])
        filter_rate2d.shape = (1,) + filter_rate2d.shape
        hyp_visc = self._vorc*filter_rate2d/self._dt
        self._hyp_visg = qgt.spec2grid(hyp_visc)
    
    def cal_barotropic(self, num_equlats=1024, rm_forcing=None, vorg=None, **args):
        if vorg is None:
            self._vorg = qgt.spec2grid(self._vorc)
        else:
            self._vorg = vorg
        if rm_forcing is not None:
            rm_forc_g = qgt.spec2grid(rm_forcing)
            mean_btpv, btpv_A, fs, As = time_mean_wave_activity_op(self._vorg, 
                num_equlats, 
                [self._hyp_visg,
                -self.params['bot_drag']*self._vorg,
                rm_forc_g], 
                beta=self.params['beta'], **args)
        else:
            mean_btpv, btpv_A, fs, As = time_mean_wave_activity_op(self._vorg, 
                num_equlats, 
                [self._hyp_visg,
                -self.params['bot_drag']*self._vorg], 
                beta=self.params['beta'], **args)
        self.lang_btpv = mean_btpv
        self.A_btpv    = btpv_A
        self.hypvis_btpv    = As[0]
        self.drag_btpv      = As[1]
        self.lang_drag_btpv = fs[1]
        if rm_forcing is not None:
            self.rm_forc = As[2]

def vorflux_1d(bt_psi):
    """
    Args:
        bt_psi: 3d numpy array with shape (time, ky, kx)
    """
    vg   = qgt.spec2grid(qgt.get_velocities(bt_psi)[1])
    vorg = qgt.spec2grid(qgt.get_vorticity(bt_psi))
    return np.mean(np.mean(vg*vorg, 0), -1)

def dfdy(fs, dy=None):
    """1th derivative of an array fs
    Args:
        fs: 1d numpy array.
        dy: spacing
    Returns:
        d:  1d numpy array with the same shape as fs
    """
    if dy is None:
        dy = 2*np.pi/fs.size
    d = np.zeros_like(fs)
    d[1:-1] = (fs[2:] - fs[:-2])/2/dy
    d[0]    = (fs[1] - fs[0])/dy
    d[-1]   = (fs[-1] - fs[-2])/dy
    return d