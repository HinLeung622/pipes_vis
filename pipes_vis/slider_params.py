from . import utils


# build a library of standard slider labels and extremes, side sorting and priority
sfh_priorities = {
    'burst':0,
    'constant':10,
    'exponential':20,
    'delayed':30,
    'lognormal':40,
    'dblplaw':50,
    'psb2':60
}

massformed_unit = r'$\log_{10}M_*/M_\odot$'
metallicity_unit = r'$Z_*/Z_\odot$'
slider_lib = {
    # various sfh parameters
    # =================== burst ===================
    'burst:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'burst:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
    'burst:tform':{
        'label':r'$t_{burst}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value],
        'side':'left',
        'priority':3
    },
    # =================== constant ==================
    'constant:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'constant:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
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
    'exponential:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'exponential:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
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
    'delayed:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'delayed:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
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
    'lognormal:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'lognormal:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
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
    'dblplaw:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'dblplaw:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
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
    # =================== psb2 ====================
    'psb2:massformed':{
        'label':massformed_unit,
        'lims':[8.5,12],
        'side':'left',
        'priority':1
    },
    'psb2:metallicity':{
        'label':metallicity_unit,
        'lims':[0.0,2.0],
        'side':'left',
        'priority':2
    },
    'psb2:told':{
        'label':r'$t_{old}$ (Gyr)',
        'lims':[0,utils.cosmo.age(0).value-2],
        'side':'left',
        'priority':3
    },
    'psb2:tau':{
        'label':r'$\tau$ (Gyr)',
        'lims':[0.1,10],
        'side':'left',
        'priority':4
    },
    'psb2:tburst':{
        'label':r'$t_{burst}$ (Gyr)',
        'lims':[utils.cosmo.age(0).value-5,utils.cosmo.age(0).value-0.3],
        'side':'left',
        'priority':5
    },
    'psb2:alpha':{
        'label':r'$\alpha$',
        'lims':[0.01,500],
        'side':'left',
        'priority':6
    },
    'psb2:beta':{
        'label':r'$\beta$',
        'lims':[0.01,500],
        'side':'left',
        'priority':7
    },
    'psb2:fburst':{
        'label':r'$f_{burst}$',
        'lims':[0.0,1.0],
        'side':'left',
        'priority':8
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

for key in slider_lib.keys():
    if ':' in key:
        for sfh_key in sfh_priorities.keys():
            if key.split(':')[0] == sfh_key:
                slider_lib[key]['priority'] += sfh_priorities[sfh_key]
