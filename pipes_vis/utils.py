import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.ndimage import median_filter

import bagpipes as pipes
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
    
    def psb_wild2020(vis_comp):
        pipes_comp = {}
        pipes_comp["age"] = cosmo.age(0).value - vis_comp['told']
        pipes_comp["tau"] = vis_comp['tau']
        pipes_comp["burstage"] = cosmo.age(0).value - vis_comp['tburst']
        pipes_comp["alpha"] = vis_comp['alpha']
        pipes_comp["beta"] = vis_comp['beta']
        pipes_comp["fburst"] = vis_comp['fburst']
        return pipes_comp
    
    def psb_twin(vis_comp):
        pipes_comp = {}
        pipes_comp["age"] = cosmo.age(0).value - vis_comp['told']
        pipes_comp["alpha1"] = vis_comp['alpha1']
        pipes_comp["beta1"] = vis_comp['beta1']
        pipes_comp["burstage"] = cosmo.age(0).value - vis_comp['tburst']
        pipes_comp["alpha2"] = vis_comp['alpha2']
        pipes_comp["beta2"] = vis_comp['beta2']
        pipes_comp["fburst"] = vis_comp['fburst']
        return pipes_comp
    
class metallicity_translate:
    def delta(vis_comp, z):
        pipes_comp = {}
        pipes_comp["metallicity_type"] = 'delta'
        pipes_comp["metallicity"] = vis_comp['metallicity']
        return pipes_comp
    
    def psb_two_step(vis_comp, z):
        pipes_comp = {}
        pipes_comp["metallicity_type"] = "two_step"
        pipes_comp["metallicity_old"] = vis_comp["metallicity_old"]
        pipes_comp["metallicity_new"] = vis_comp["metallicity_burst"]
        pipes_comp["metallicity_step_age"] = cosmo.age(z).value - vis_comp['tburst']
        return pipes_comp
    
    def two_step(vis_comp, z):
        pipes_comp = {key:vis_comp[key] for key in ["metallicity_type", "metallicity_old", "metallicity_new"]}
        pipes_comp["metallicity_step_age"] = cosmo.age(z).value - vis_comp['metallicity_tstep']
        return pipes_comp
    
    def psb_linear_step(vis_comp, z):
        pipes_comp = {}
        pipes_comp["metallicity_type"] = "linear_step"
        pipes_comp["metallicity_zero"] = vis_comp["metallicity_zero"] + \
            vis_comp["metallicity_slope"]*(cosmo.age(0).value-cosmo.age(z).value)
        pipes_comp["metallicity_slope"] = vis_comp["metallicity_slope"]
        pipes_comp["metallicity_new"] = vis_comp["metallicity_burst"]
        pipes_comp["metallicity_step_age"] = cosmo.age(z).value - vis_comp['tburst']
        return pipes_comp
    
    def linear_step(vis_comp, z):
        pipes_comp = {key:vis_comp[key] for key in ["metallicity_type", "metallicity_slope", "metallicity_new"]}
        pipes_comp["metallicity_zero"] = vis_comp["metallicity_zero"] + \
            vis_comp["metallicity_slope"]*(cosmo.age(0).value-cosmo.age(z).value)
        pipes_comp["metallicity_step_age"] = cosmo.age(z).value - vis_comp['metallicity_tstep']
        return pipes_comp

class get_ceh_array:
    """
    Evaluates the metallicity values at a list of ages (in normal age) given the 
    metallicity model choice and model parameters.
    """
    def delta(ages, sfh_dict):
        return np.ones(len(ages))*sfh_dict['metallicity']
    
    def two_step(ages, sfh_dict):
        pre_step_ind = np.where(ages < sfh_dict['metallicity_tstep'])
        post_step_ind = np.isin(np.arange(len(ages)), pre_step_ind, invert=True)
        ceh = np.zeros(len(ages))
        ceh[pre_step_ind] = sfh_dict['metallicity_old']
        ceh[post_step_ind] = sfh_dict['metallicity_new']
        return ceh
    
    def psb_two_step(ages, sfh_dict):
        pre_step_ind = np.where(ages < sfh_dict['tburst'])
        post_step_ind = np.isin(np.arange(len(ages)), pre_step_ind, invert=True)
        ceh = np.zeros(len(ages))
        ceh[pre_step_ind] = sfh_dict['metallicity_old']
        ceh[post_step_ind] = sfh_dict['metallicity_burst']
        return ceh
    
    def linear_step(ages, sfh_dict):
        pre_step_ind = np.where(ages < sfh_dict['metallicity_tstep'])
        post_step_ind = np.isin(np.arange(len(ages)), pre_step_ind, invert=True)
        ceh = np.zeros(len(ages))
        ceh[pre_step_ind] = sfh_dict['metallicity_zero'] + \
            sfh_dict['metallicity_slope']*(cosmo.age(0).value-ages[pre_step_ind])
        ceh[post_step_ind] = sfh_dict['metallicity_new']
        return ceh
    
    def psb_linear_step(ages, sfh_dict):
        pre_step_ind = np.where(ages < sfh_dict['tburst'])
        post_step_ind = np.isin(np.arange(len(ages)), pre_step_ind, invert=True)
        ceh = np.zeros(len(ages))
        ceh[pre_step_ind] = sfh_dict['metallicity_zero'] + \
            sfh_dict['metallicity_slope']*(cosmo.age(0).value-ages[pre_step_ind])
        ceh[post_step_ind] = sfh_dict['metallicity_burst']
        return ceh

def make_sfh_comp(sfh_type, z0_comp):
    """Creates a sfh dictionary at 0 redshift, but in lookback time"""
    sfh_comp = getattr(sfh_translate, sfh_type)(z0_comp)
    sfh_comp['massformed'] = z0_comp['massformed']
    if "metallicity_type" in z0_comp.keys():
        metallicity_comp = getattr(metallicity_translate, z0_comp["metallicity_type"])(z0_comp, 0)
    else:
        metallicity_comp = getattr(metallicity_translate, 'delta')(z0_comp, 0)
    sfh_comp.update(metallicity_comp)
    return sfh_comp

def parametric_to_custom_sfh(sfh_type, sfh_dict, ages):
    """
    Using the Bagpipes models.sfh object, to extract an array of age (lb time)
    vs sfr. Then interp this SFH onto a given age list.
    """
    model_components = {}
    model_components[sfh_type] = make_sfh_comp(sfh_type, sfh_dict)
    model_components['redshift'] = 0

    # pass the temp input dictionary to the SFH module in bagpipes to build an SFH
    pipes_sfh = pipes.models.star_formation_history(model_components)
    to_interp_ages = cosmo.age(0).value*10**9 - pipes_sfh.ages[::-1]
    sfh = pipes_sfh.sfh[::-1]
    
    return np.interp(ages, to_interp_ages/10**9, sfh)

def create_sfh(params):
    """
    Using the visualizer input dictionary to create custom SFH arrays for each
    SFH component fixed at z=0 with Bagpipes, then translate that array to the 
    given z, and return 3 things: 1. total log10M before given z, 2. SFH array 
    used to plot SFH and 3. a dictionary of custom SFH dictionaries that is to 
    be poured into bagpipes' model components
    """
    #print('')
    pipes_sfh_funcs = dir(pipes.models.star_formation_history)
    #for key in params.keys():
    #    if key in pipes_sfh_funcs or key[:-1] in pipes_sfh_funcs:
    #        print(params[key])
    
    z = params["redshift"]
    total_M_before_z = 0
    total_sfh = None
    pipes_sfh_dict = {}
    # create a temporary bagpipes input dictionary that only has the SFH-related elements, z forced to be 0
    pipes_sfh_funcs = dir(pipes.models.star_formation_history)
    custom_count = 0
    for key in params.keys():
        if key in pipes_sfh_funcs or key[:-1] in pipes_sfh_funcs:
            custom_count += 1
            # bagpipes is only able to deal with 9 sfh components of the same type, so error is raised when
            # more than 9 sfh components are set (they all have to use the custom component to pass to bagpipes)
            if custom_count > 9:
                raise NameError("Too many initiated SFH components, exceeds bagpipes' limits(9).")
            model_components = {}
            if key in pipes_sfh_funcs:
                model_components[key] = make_sfh_comp(key, params[key])
            elif key[:-1] in pipes_sfh_funcs:
                model_components[key] = make_sfh_comp(key[:-1], params[key])
            model_components['redshift'] = 0
    
            # pass the temp input dictionary to the SFH module in bagpipes to build an SFH
            pipes_sfh = pipes.models.star_formation_history(model_components)
            lb_time = pipes_sfh.ages
            sfh = pipes_sfh.sfh
            
            # insert an interpolated datapoint in bagpipes generated SFH array at redshift z, 
            # to prevent jittering of model spectrum due to such interpolate issues
            lb_time_of_z = (cosmo.age(0).value - cosmo.age(z).value)*10**9
            sfr_at_z = np.interp(lb_time_of_z, lb_time, sfh)
            
            # calculate the mass formed up to the given redshift z
            age = cosmo.age(0).value*10**9 - lb_time[::-1]
            sfh = sfh[::-1]
            age_before_z = age[age<cosmo.age(z).value*10**9]
            sfh_before_z = sfh[age<cosmo.age(z).value*10**9]
            M_before_z = np.trapz(y=sfh_before_z, x=age_before_z)
            
            #create custom sfh array that is before z for input to model components
            #requirements: column 0 = lookback time at z in years, column 1 = SFR in Msuns per year
            lb_time_before_z = cosmo.age(z).value*10**9 - age_before_z[::-1]
            lb_sfh_before_z = sfh_before_z[::-1]
            lb_time_before_z = np.insert(lb_time_before_z, 0, 0)        # the insert
            lb_sfh_before_z = np.insert(lb_sfh_before_z, 0, sfr_at_z)   # the insert
            custom_sfh = np.array([lb_time_before_z, lb_sfh_before_z]).T
            
            # now add this sfh component's sfh to the total
            if total_sfh is not None:
                total_sfh += sfh
            else:
                total_sfh = sfh
            # build up the total mass formed from the different SFH components
            total_M_before_z += M_before_z
            # put everything into a dictionary that will feed directly into bagpipes
            custom = {}
            custom['massformed'] = np.log10(M_before_z)
            # check metallicity and translate
            if "metallicity_type" in params[key].keys():
                metallicity_comp = getattr(metallicity_translate, params[key]["metallicity_type"])(params[key], z)
            else:
                metallicity_comp = getattr(metallicity_translate, 'delta')(params[key], z)
            custom.update(metallicity_comp)
            
            #print(custom)
            
            # set the sfh array to non zero throughout to prevent bagpipes going crazy
            if M_before_z == 0:
                custom_sfh[:,1] = np.ones(len(custom_sfh))
            custom['history'] = custom_sfh
            pipes_sfh_dict['custom'+str(custom_count)] = custom
    
    return np.log10(total_M_before_z), np.array([age/10**9, total_sfh]), pipes_sfh_dict

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
    
def make_pipes_components(params, sfh_dict):
    """ Create a Bagpipes component dictionary (custom SFHs) from visualizer inputs """
    #print(params)
    #print(sfh_dict)
    model_components = sfh_dict.copy()
    model_components["redshift"] = params["redshift"]
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

def shift_index(ind_dict, redshift):
    """ Shifts indices' wavelengths from rest frame to observed frame """
    out_dict = {}
    for key in ind_dict:
        if key == 'feature':
            out_dict['feature'] = np.array(ind_dict['feature'])*(1+redshift)
        elif key == 'continuum':
            out_dict['continuum'] = []
            for i in range(len(ind_dict['continuum'])):
                out_dict['continuum'].append(np.array(ind_dict['continuum'][i])*(1+redshift))
        else:
            out_dict[key] = ind_dict[key]
    return out_dict

def get_ceh(ages, params, model):
    """
    Get the mass-weighted metallicity array at a specified list of ages (lb time)
    for one or more SFH component
    """
    sfh_comp_keys = []
    pipes_sfh_funcs = dir(pipes.models.star_formation_history)
    for key in params.keys():
        if key in pipes_sfh_funcs or key[:-1] in pipes_sfh_funcs:
            sfh_comp_keys.append(key)
            
    if len(sfh_comp_keys) == 1:
        if 'metallicity_type' in params[sfh_comp_keys[0]].keys():
            zmet_evo = getattr(get_ceh_array, params[sfh_comp_keys[0]]['metallicity_type'])(
                ages, params[sfh_comp_keys[0]])
        else:
            zmet_evo = getattr(get_ceh_array, 'delta')(ages, params[sfh_comp_keys[0]])
            
        zmet_evo[zmet_evo<0] = 0
        return zmet_evo
            
    else:
        zmet_evo = np.zeros([len(sfh_comp_keys), len(ages)])
        sfh_array = np.zeros(np.shape(zmet_evo))
        for i,key in enumerate(sfh_comp_keys):
            if 'metallicity_type' in params[key].keys():
                zmet_evo[i,:] = getattr(get_ceh_array, params[key]['metallicity_type'])(ages, params[key])
            else:
                zmet_evo[i,:] = getattr(get_ceh_array, 'delta')(ages, params[key])
                
            if key in pipes_sfh_funcs:
                sfh_array[i,:] = parametric_to_custom_sfh(key, params[key], ages)
            elif key[:-1] in pipes_sfh_funcs:
                sfh_array[i,:] = parametric_to_custom_sfh(key[:-1], params[key], ages)
            zmet_evo[i,:][zmet_evo[i,:]<0] = 0
        
        sfh_sum = np.sum(sfh_array, axis=0)
        zmet_evo_1d = np.sum(zmet_evo*sfh_array, axis=0) / sfh_sum
        # mask all 0 SFR bins into np.nan to prevent plotting anomalies
        zmet_evo_1d[sfh_sum<=0] = np.nan
        return zmet_evo_1d

cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
