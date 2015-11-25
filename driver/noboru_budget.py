filename_prefix = 'Nov4_Sc2.0_drag5e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg110' %filename_prefix
psic = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'psi'))

params = {'beta': 791.571747206000,
          'F'   : 395.785873603000,
          'bot_drag' : 1.989436788650000E-002}

finite_wave_diag = TwoLayerWaveActivity(psic, params)
finite_wave_diag.cal_upperlayer()
finite_wave_diag.cal_lowerlayer()