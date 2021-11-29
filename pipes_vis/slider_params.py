from . import utils


# build a library of standard slider labels and extremes, side sorting and priority
sfh_priorities = {
    'burst':0,
    'constant':10,
    'exponential':20,
    'delayed':30,
    'lognormal':40,
    'dblplaw':50,
    'psb_wild2020':60,
    'psb_twin':70
}

# commonly used among all SFHs
massformed_dict = {
    'label':r'$\log_{10}M_*/M_\odot$',
    'lims':[8.5,12],
    'side':'left',
    'priority':1
    }

# library for all possible metallicity parameters, added to each possible SFH component
metallicity_lib = {
    # constant metallicity
    'metallicity':{
        'label':r'$Z_*/Z_\odot$',
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2.0
        },
    # psb_two_step and two_step
    'metallicity_old':{
        'label':r'$Z_{old}/Z_\odot$',
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2.1
        },
    # psb_two_step, and psb_linear_step
    'metallicity_burst':{
        'label':r'$Z_{burst}/Z_\odot$',
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2.2
        },
    # two_step
    'metallicity_new':{
        'label':r'$Z_{new}/Z_\odot$',
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2.3
        },
    # psb_linear_step
    'metallicity_slope':{
        'label':r'$Z_*/Z_\odot \; Gyr^{-1}$',
        'lims':[-1.0,1.0],
        'side':'left',
        'priority':2.4
        },
    # psb_linear_step
    'metallicity_zero':{
        'label':r'$Z_0/Z_\odot$',
        'lims':[0.0,5.0],
        'side':'left',
        'priority':2.5
        },
    # two_step
    'metallicity_tstep':{
        'label':r'$t_{Zstep}$ (Gyr)',
        'lims':[0.0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':2.6
        }
    }

slider_lib = {
    # various sfh parameters
    # =================== burst ===================
    'burst:tform':{
        'label':r'$t_{burst}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':3
    },
    # =================== constant ==================
    'constant:tform':{
        'label':r'$t_{form}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':3
    },
    'constant:tend':{
        'label':r'$t_{end}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':4
    },
    # =================== exponential ====================
    'exponential:tform':{
        'label':r'$t_{form}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value-1],
        'side':'left',
        'priority':3
    },
    'exponential:tau':{
        'label':r'$\tau$ (Gyr)',
        'lims':[0.1,10],
        'side':'left',
        'priority':4
    },
    # =================== delayed ====================
    'delayed:tform':{
        'label':r'$t_{form}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value-1],
        'side':'left',
        'priority':3
    },
    'delayed:tau':{
        'label':r'$\tau$ (Gyr)',
        'lims':[0.1,10],
        'side':'left',
        'priority':4
    },
    # =================== lognormal ====================
    'lognormal:tmax':{
        'label':r'$t_{max}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':3
    },
    'lognormal:fwhm':{
        'label':'FWHM (Gyr)',
        'lims':[0.1,10],
        'side':'left',
        'priority':4
    },
    # =================== dblplaw ====================
    'dblplaw:tau':{
        'label':r'$\tau\;/\;t_{max}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':3
    },
    'dblplaw:alpha':{
        'label':r'$\alpha$',
        'lims':[0.01,500],
        'side':'left',
        'priority':4
    },
    'dblplaw:beta':{
        'label':r'$\beta$',
        'lims':[0.01,500],
        'side':'left',
        'priority':5
    },
    # =================== psb_wild2020 ====================
    'psb_wild2020:told':{
        'label':r'$t_{old}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value-2],
        'side':'left',
        'priority':3
    },
    'psb_wild2020:tau':{
        'label':r'$\tau$ (Gyr)',
        'lims':[0.1,10],
        'side':'left',
        'priority':4
    },
    'psb_wild2020:tburst':{
        'label':r'$t_{burst}$ (Gyr)',
        'lims':[utils.cosmo.age(0).value-5,utils.cosmo.age(0).value-0.3],
        'side':'left',
        'priority':5
    },
    'psb_wild2020:alpha':{
        'label':r'$\alpha$',
        'lims':[0.01,500],
        'side':'left',
        'priority':6
    },
    'psb_wild2020:beta':{
        'label':r'$\beta$',
        'lims':[0.01,500],
        'side':'left',
        'priority':7
    },
    'psb_wild2020:fburst':{
        'label':r'$f_{burst}$',
        'lims':[0.0,1.0],
        'side':'left',
        'priority':8
    },
    # =================== psb_twin ===================
    'psb_twin:told':{
        'label':r'$t_{old}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value-2],
        'side':'left',
        'priority':3
    },
    'psb_twin:alpha1':{
        'label':r'$\alpha_1$',
        'lims':[0.01,5],
        'side':'left',
        'priority':4
    },
    'psb_twin:beta1':{
        'label':r'$\beta_1$',
        'lims':[0.01,500],
        'side':'left',
        'priority':5
    },
    'psb_twin:tburst':{
        'label':r'$t_{burst}$ (Gyr)',
        'lims':[utils.cosmo.age(0).value-5,utils.cosmo.age(0).value-0.3],
        'side':'left',
        'priority':6
    },
    'psb_twin:alpha2':{
        'label':r'$\alpha_2$',
        'lims':[0.01,500],
        'side':'left',
        'priority':7
    },
    'psb_twin:beta2':{
        'label':r'$\beta_2$',
        'lims':[0.01,500],
        'side':'left',
        'priority':8
    },
    'psb_twin:fburst':{
        'label':r'$f_{burst}$',
        'lims':[0.0,1.0],
        'side':'left',
        'priority':9
    },
    # =================== now right hand sides ===================
    'redshift':{
        'label':'redshift',
        'lims':[0,2],
        'side':'right',
        'priority':1
    },
    'dust:eta':{
        'label':r'$\eta_{dust}$',
        'lims':[1,4],
        'side':'right',
        'priority':2
    },
    'dust:Av':{
        'label':r'$A_V$',
        'lims':[0,2],
        'side':'right',
        'priority':3
    },
    'dust:n':{
        'label':r'$n_{CF00}$',
        'lims':[0.3,2.7],
        'side':'right',
        'priority':4
    },
    'dust:delta':{
        'label':r'$\delta_{Salim}$',
        'lims':[-1,1],
        'side':'right',
        'priority':5
    },
    'dust:B':{
        'label':r'$B_{Salim}$',
        'lims':[0,5],
        'side':'right',
        'priority':6
    },
    'dust:qpah':{
        'label':r'$q_{PAH}$ (%)',
        'lims':[0,5],
        'side':'right',
        'priority':7
    },
    'dust:umin':{
        'label':r'$U_{min}$',
        'lims':[0,10],
        'side':'right',
        'priority':8
    },
    'dust:gamma':{
        'label':r'$\gamma$',
        'lims':[0,0.1],
        'side':'right',
        'priority':9
    },
    'nebular:logU':{
        'label':r'$\log_{10}U$',
        'lims':[-10,-2],
        'side':'right',
        'priority':10
    },
    't_bc':{
        'label':r'$t_{bc}$ (Gyr)',
        'lims':[0.001,0.1],
        'side':'right',
        'priority':11
    },
    'veldisp':{
        'label':r'veldisp$\sigma$ (km/s)',
        'lims':[0,400],
        'side':'right',
        'priority':12
    }
}

for sfh_key in sfh_priorities.keys():
    slider_lib[sfh_key+':massformed'] = massformed_dict.copy()
    for metallicity_key in metallicity_lib:
        slider_lib[sfh_key+':'+metallicity_key] = metallicity_lib[metallicity_key].copy()

for key in slider_lib.keys():
    if ':' in key:
        for sfh_key in sfh_priorities.keys():
            if key.split(':')[0] == sfh_key:
                slider_lib[key]['priority'] += sfh_priorities[sfh_key]
