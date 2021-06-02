import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.ndimage import median_filter

try:
    import bagpipes as pipes
except ImportError:
    print('BAGPIPES not installed')
from . import override_config
override_config.override_config(pipes)


class sfh_translate:
    def burst(vis_comp):
        pipes_comp = {}
        pipes_comp['age'] = cosmo.age(0).value - vis_comp['tform']
        return pipes_comp
    
    def constant(vis_comp):
        pipes_comp = {}
        pipes_comp["age_max"] = cosmo.age(0).value - vis_comp['tform']
        pipes_comp["age_min"] = cosmo.age(0).value - vis_comp['tend']
        return pipes_comp
    
    def exponential(vis_comp):
        pipes_comp = {}
        pipes_comp["age"] = cosmo.age(0).value - vis_comp['tform']
        pipes_comp["tau"] = vis_comp['tau']
        return pipes_comp
    
    def delayed(vis_comp):
        pipes_comp = {}
        pipes_comp["age"] = cosmo.age(0).value - vis_comp['tform']
        pipes_comp["tau"] = vis_comp['tau']
        return pipes_comp
    
    def lognormal(vis_comp):
        pipes_comp = {}
        pipes_comp["tmax"] = vis_comp['tmax']
        pipes_comp["fwhm"] = vis_comp['fwhm']
        return pipes_comp
    
    def dblplaw(vis_comp):
        pipes_comp = {}
        pipes_comp["tau"] = vis_comp['tau']
        pipes_comp["alpha"] = vis_comp['alpha']
        pipes_comp["beta"] = vis_comp['beta']
        return pipes_comp
    
    def psb2(vis_comp):
        pipes_comp = {}
        pipes_comp["age"] = cosmo.age(0).value - vis_comp['told']
        pipes_comp["tau"] = vis_comp['tau']
        pipes_comp["burstage"] = cosmo.age(0).value - vis_comp['tburst']
        pipes_comp["alpha"] = vis_comp['alpha']
        pipes_comp["beta"] = vis_comp['beta']
        pipes_comp["fburst"] = vis_comp['fburst']
        return pipes_comp

def create_sfh(params):
    """
    Using the visualizer input dictionary to create a custom SFH array fixed at 
    z=0 with Bagpipes, then translate that array to the given z, and return 3
    things: 1. total log10M before given z, 2. SFH array used to plot SFH and 
    3. custom SFH that is to be poured into bagpipes' custom sfh module
    """
    z = params["redshift"]
    # create a temporary bagpipes input dictionary that only has the SFH-related elements, z forced to be 0
    model_components = {}
    for key in params.keys():
        if key in dir(pipes.models.star_formation_history):
            sfh_comp = getattr(sfh_translate, key)(params[key])
            sfh_comp['massformed'] = params[key]['massformed']
            sfh_comp['metallicity'] = params['metallicity']
            model_components[key] = sfh_comp
        elif key[:-1] in dir(pipes.models.star_formation_history):
            sfh_comp = getattr(sfh_translate, key[:-1])(params[key])
            sfh_comp['massformed'] = params[key]['massformed']
            sfh_comp['metallicity'] = params['metallicity']
            model_components[key] = sfh_comp
    model_components['redshift'] = 0
    
    # pass the temp input dictionary to the SFH module in bagpipes to build an SFH
    pipes_sfh = pipes.models.star_formation_history(model_components)
    lb_time = pipes_sfh.ages
    sfh = pipes_sfh.sfh
    
    # insert an interpolated datapoint in bagpipes generated SFH array at redshift z, 
    # to prevent jittering of model spectrum due to such interpolate issues
    lb_time_of_z = (cosmo.age(0).value - cosmo.age(z).value)*10**9
    sfr_at_z = np.interp(lb_time_of_z, lb_time, sfh)
    
    #calculate the mass formed up to the given redshift z
    age = cosmo.age(0).value*10**9 - lb_time[::-1]
    sfh = sfh[::-1]
    age_before_z = age[age<cosmo.age(z).value*10**9]
    sfh_before_z = sfh[age<cosmo.age(z).value*10**9]
    logM_before_z = np.log10(np.trapz(y=sfh_before_z, x=age_before_z))
    
    #create custom sfh array that is before z for input to model components
    #requirements: column 0 = lookback time at z in years, column 1 = SFR in Msuns per year
    lb_time_before_z = cosmo.age(z).value*10**9 - age_before_z[::-1]
    lb_sfh_before_z = sfh_before_z[::-1]
    lb_time_before_z = np.insert(lb_time_before_z, 0, 0)        # the insert
    lb_sfh_before_z = np.insert(lb_sfh_before_z, 0, sfr_at_z)   # the insert
    custom_sfh = np.array([lb_time_before_z, lb_sfh_before_z]).T
    
    return logM_before_z, np.array([age/10**9, sfh]), custom_sfh

'''
def running_median(x,y, width=150, indexing=None):
    """
    Calculates a rolling/running median for spectrum.
    The search for wavelengths within the width will not be performed if pre-
    calculated indexing is provided.
    """
    # width = width of wavelength considered in each median taken
    medians = np.zeros(len(y))
    if indexing is not None:
        for i,yi in enumerate(y):
            y_sec = y[indexing[i]]
            medians[i] = np.median(y_sec)
    else:
        for i,yi in enumerate(y):
            y_sec = y[np.where((x >= x[i]-width/2) & (x <= x[i]+width/2))]
            medians[i] = np.median(y_sec)
    return medians
'''

def running_median(x,y, width=150):
    """ Calculates a rolling/running median for spectrum """
    # width = width of wavelength considered in each median taken
    x_gaps = np.diff(x)
    # get indices of transitions where the spacing between concecutive x entries change
    gap_diffs = np.diff(x_gaps)
    gap_change_ind = np.where(gap_diffs>0.01)[0]+2
    gap_change_ind = np.insert(gap_change_ind, len(gap_change_ind), len(y))
    #print(gap_change_ind, len(gap_change_ind))
    medians = np.zeros(len(y))
    for i,indi in enumerate(gap_change_ind):
        this_gap = x[indi-1] - x[indi-2]
        window_size = int(np.floor(width/this_gap))
        if i != 0:
            l_indi = gap_change_ind[i-1]
            left_ind = np.where((x >= x[l_indi]-width/2) & (x <= x[l_indi-1]))
            x_left = x[left_ind]
            y_left = y[left_ind]
            x_left_ = np.arange(1,width/2/this_gap+1)*this_gap
            x_left_ = x[l_indi] - x_left_[::-1]
            y_left_ = np.interp(x_left_, x_left, y_left)
            
            if i == len(gap_change_ind)-1:
                # last bin
                y_sec = np.concatenate([y_left_, y[l_indi:]])
                #x_sec = np.concatenate([x_left_, x[l_indi:]])
                #print('last bin')
                #print('main body', x[l_indi:][0], x[l_indi:][-1])
                #print('x_left_', x_left_[0], x_left_[-1])
                #print(x_sec[0], x_sec[-1], set(np.diff(x_sec)))
                median_sec = median_filter(y_sec, size=window_size)[len(x_left_):]
                medians[l_indi:] += median_sec
                
            else:
                # middle bins
                right_ind = np.where((x >= x[indi-1]) & (x <= x[indi]+width/2))
                x_right = x[right_ind]
                y_right = y[right_ind]
                x_right_ = np.arange(1,width/2/this_gap+1)*this_gap + x[indi-1]
                y_right_ = np.interp(x_right_, x_right, y_right)
                y_sec = np.concatenate([y_left_, y[l_indi:indi], y_right_])
                #x_sec = np.concatenate([x_left_, x[l_indi:indi], x_right_])
                #print(f'{i}th bin')
                #print('main body', x[l_indi:indi][0], x[l_indi:indi][-1])
                #print('x_left_', x_left_[0], x_left_[-1])
                #print('x_right_', x_right_[0], x_right[-1])
                #print(x_sec[0], x_sec[-1], set(np.diff(x_sec)))
                median_sec = median_filter(y_sec, size=window_size)[len(x_left_):-len(x_right_)]
                medians[l_indi:indi] += median_sec
                
        else:
            # first bin
            right_ind = np.where((x >= x[indi-1]) & (x <= x[indi]+width/2))
            x_right = x[right_ind]
            y_right = y[right_ind]
            x_right_ = np.arange(1,width/2/this_gap+1)*this_gap + x[indi-1]
            y_right_ = np.interp(x_right_, x_right, y_right)
            y_sec = np.concatenate([y[:indi], y_right_])
            #x_sec = np.concatenate([x[:indi], x_right_])
            #print('first bin')
            #print('main body', x[:indi][0], x[:indi][-1])
            #print('x_right_', x_right_[0], x_right_[-1])
            #print(x_sec[0], x_sec[-1], set(np.diff(x_sec)))
            median_sec = median_filter(y_sec, size=window_size)[:-len(x_right_)]
            medians[:indi] += median_sec
            
    return medians
    
def make_pipes_components(params, input_logM, sfh):
    """ Create a Bagpipes component dictionary (custom SFH) from visualizer inputs """
    custom = {}
    custom["massformed"] = input_logM
    custom["metallicity"] = params["metallicity"]
    custom["history"] = sfh

    model_components = {}
    model_components["redshift"] = params["redshift"]
    model_components["custom"] = custom
    if "dust" in params.keys():
        model_components["dust"] = params["dust"]
    if "nebular" in params.keys():
        model_components["nebular"] = params["nebular"]
    if "t_bc" in params.keys():
        model_components["t_bc"] = params["t_bc"]
    if "veldisp" in params.keys():
        if params["veldisp"]>0:
            model_components["veldisp"] = params["veldisp"]
    return model_components

cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
