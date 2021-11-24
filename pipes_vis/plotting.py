import numpy as np

import bagpipes as pipes
from . import override_config
override_config.override_config(pipes)

from . import utils


def add_bp_sfh(sfh, ax, color="black"):
    """ Adds a SFH black line in the given axis, in Bagpipes style """
    zvals=[0, 0.5, 1, 2, 4, 10]
    sfh_line, = ax.plot(sfh[0], sfh[1],
            color=color, zorder=4, lw=2, alpha=1, ls="-", label=None)

    # Set limits
    ax.set_xlim([utils.cosmo.age(0).value, 0.])
    ax.set_ylim(bottom=0.)

    # Add redshift axis along the top
    pipes.plotting.add_z_axis(ax, zvals=zvals)

    # Add labels
    ax.set_ylabel("$\\mathrm{SFR\\ /\\ M_\\odot\\ \\mathrm{yr}^{-1}}$")
    ax.set_xlabel("$\\mathrm{Age\\ of\\ Universe\\ /\\ Gyr}$")
    return sfh_line

def add_sfh_plot(params, sfh, input_logM, ax, sfh_color="black", sfh_label=None, 
                 z_line_color='red', z_line_label='obs redshift'):
    """
    Adds the elements in the sfh panel to the given axis,
    returns plot elements that can be later updated by sliders
    """
    z = params["redshift"]
    sfh_line = add_bp_sfh(sfh, ax, color=sfh_color)
    age_at_z = utils.cosmo.age(z).value
    z_line, = ax.plot([age_at_z,age_at_z], [0,ax.get_ylim()[1]*1000], color=z_line_color, 
                         label=z_line_label)
    z_text = ax.annotate(str(np.round(z,3)), [age_at_z, 0.92*ax.get_ylim()[1]], color=z_line_color)
    ax.legend(loc='upper right')
    input_logM_text = ax.text(0.985,0.88,r'pre-obs $\log_{10}M_*$ = '+str(np.round(input_logM,2)),
                   transform=ax.transAxes,fontsize=14, horizontalalignment='right',
                   verticalalignment='top', zorder=10)
    
    bad_sfh_text = ax.text(0.5, 0.5, 'Oops! Bad SFH inputs, no SFH displayed', color='red',
                           fontsize=24, alpha=0, horizontalalignment='center', 
                           verticalalignment='center', transform=ax.transAxes, zorder=100)
    
    return sfh_line, z_line, z_text, input_logM_text, bad_sfh_text

def add_bp_spectrum(spectrum, ax, zorder=4, z_non_zero=True, color="sandybrown",
                    ymax=None, lw=2., label=None, alpha=1, sub=False):
    """ 
    Add a spectrum to the passed axes. Adds errors if they are
    included in the spectrum object as a third column. 
    """
    # Sort out axis limits
    if not ymax:
        ymax = 1.05*np.max(spectrum[:, 1])

    y_scale = int(np.log10(ymax))-1

    ax.set_ylim(0., ymax*10**-y_scale)
    ax.set_xlim(spectrum[0, 0], spectrum[-1, 0])

    # Plot the data
    spec_line, = ax.plot(spectrum[:, 0], spectrum[:, 1]*10**-y_scale,
            color=color, zorder=zorder, lw=lw, label=label, alpha=alpha)

    # Sort out x tick locations
    pipes.plotting.auto_x_ticks(ax)
    if sub == True:
        ax.set_xscale('log')
        ax.set_yticks([])

    # Sort out axis labels.
    if sub == False:
        pipes.plotting.auto_axis_label(ax, y_scale, z_non_zero=z_non_zero)

    return y_scale, spec_line

def add_main_spec(model, ax, spec_lim, median_width=150, color="sandybrown", 
                  continuum_color='black', continuum_label='continuum'):
    """
    Adds the elements in the main spectrum plot to the given axis,
    returns plot elements that can be later updated by sliders and y scale of spectrum
    """
    run_med = utils.running_median(model.spectrum[:,0], model.spectrum[:,1], width=median_width)
    zoom_in_spec = model.spectrum[np.where((model.spectrum[:,0] >= spec_lim[0]) & 
                                           (model.spectrum[:,0] <= spec_lim[1]))]
    y_scale_spec,spec_line = add_bp_spectrum(zoom_in_spec, ax, color=color)
    run_med_line, = ax.plot(model.spectrum[:, 0], run_med*10**-y_scale_spec, color='k', alpha=0.3, zorder=5,
                            label='continuum')
    ax.legend(loc='upper right')
    overflow_text = ax.text(0.5, 0.5, 'Fluxes too low, y-axis scale inaccurate!!', color='red',
                            fontsize=24, alpha=0, horizontalalignment='center', 
                            verticalalignment='center', transform=ax.transAxes, zorder=100)
    
    return spec_line, run_med_line, overflow_text, y_scale_spec

def update_spec(spectrum, ax, spec_line, sub=False, overflow_text=None, change_xlims=False):
    """ updates a given spectrum line plot from a given ax with a given new spectrum """
    if overflow_text is not None:
        overflow_text.set_alpha(0.0)
    spec_ymax = 1.05*np.max(spectrum[:, 1])
    try:
        spec_y_scale = int(np.log10(spec_ymax))-1
        spec_line.set_ydata(spectrum[:, 1]*10**-spec_y_scale)
        if change_xlims == True:
            spec_line.set_xdata(spectrum[:, 0])
            ax.set_ylim(0., spec_ymax*10**-spec_y_scale)
            ax.set_xlim(spectrum[0, 0], spectrum[-1, 0])
            pipes.plotting.auto_x_ticks(ax)
        ax.set_ylim(0., spec_ymax*10**-spec_y_scale)
        if sub == False:
            pipes.plotting.auto_axis_label(ax, spec_y_scale, z_non_zero=True)
        return spec_y_scale
    except OverflowError:
        if overflow_text is not None:
            overflow_text.set_alpha(1.0)

def add_residual(model, ax, spec_lim, median_width=150, color="sandybrown"):
    """
    Adds the elements in the residual plot to the given axis,
    returns plot elements that can be later updated by sliders and y scale of spectrum
    """
    run_med = utils.running_median(model.spectrum[:,0], model.spectrum[:,1], width=median_width)
    residual = model.spectrum[:,1] / run_med
    res_line, = ax.plot(model.spectrum[:, 0], residual, color=color, lw=2, zorder=1)
    # value guidelines and labels
    ax.axhline(1, color="black", ls="--", lw=1, zorder=0)
    ax.axhline(1.5, color="black", ls=":", lw=1, zorder=0)
    ax.annotate('1.5x', [0.98*(spec_lim[1]-spec_lim[0])+spec_lim[0], 1.5], ha='center', va='center')
    ax.axhline(0.5, color="black", ls=":", lw=1, zorder=0)
    ax.annotate('0.5x', [0.98*(spec_lim[1]-spec_lim[0])+spec_lim[0], 0.5], ha='center', va='center')
    ax.set_xlim(spec_lim)
    pipes.plotting.auto_x_ticks(ax)
    pipes.plotting.auto_axis_label(ax, -1, z_non_zero=True)
    ax.set_ylabel('flux/\ncontinuum')
    # rescale the y axis to be determined only by the residuals in frame
    in_range_res = residual[np.where((model.spectrum[:,0] >= spec_lim[0]) & 
                                     (model.spectrum[:,0] <= spec_lim[1]))]
    res_span = max(in_range_res) - min(in_range_res)
    ax.set_ylim([min(in_range_res)-0.1*res_span, max(in_range_res)+0.1*res_span])
    
    return res_line

def add_index_spectrum(ind_dict, spectrum, ind_val, redshift, y_scale, 
                       color_continuum='lightgray', color_feature='sandybrown', alpha=0.2):
    """ 
    Adds the elements in line index plots to a given axis, 
    returns plot elements that can later be updated by sliders and y scale of spectrum
    """
    index_ax = ind_dict['ax']
    ind_dict_ = utils.shift_index(ind_dict, redshift)
    if ind_dict_['type'] in ['EW', 'break']:
        continuum_flat = [j for sub in ind_dict_['continuum'] for j in sub]
        index_range = [min(continuum_flat), max(continuum_flat)]
        spec_range = np.where((spectrum[:,0] > index_range[0]) & (spectrum[:,0] < index_range[1]))
        ind_line = index_ax.plot(spectrum[spec_range][:,0], spectrum[spec_range][:,1]*10**-y_scale, color='k')
        ind_dict['line'] = ind_line[0]

        ylims = index_ax.get_ylim()
        con_polys = []
        for con in ind_dict_['continuum']:
            con_poly = index_ax.fill_between(con, [0,0], [20,20], color='lightgray', alpha=0.2)
            con_polys.append(con_poly)
        ind_dict['con_polys'] = con_polys
        if ind_dict_['type'] == 'EW':
            feature_poly = index_ax.fill_between(ind_dict_['feature'], [0,0], [20,20], color='sandybrown', alpha=0.2)
            ind_dict['feature_poly'] = feature_poly
        index_ax.set_ylim(ylims)
        index_ax.set_xlim(index_range)
        index_ax.tick_params(direction="in")
        ind_text = index_ax.text(0.95, 1.0, ind_dict_['name']+'='+str(np.round(ind_val,2)), 
                                 ha='right', va='bottom', transform=index_ax.transAxes)
    else:
        index_ax.set_xticks([])
        ind_text = index_ax.text(0.5, 0.5, ind_dict_['name']+'='+str(np.round(ind_val,2)), 
                                 ha='center', va='center', transform=index_ax.transAxes)
    
    ind_dict['text'] = ind_text
    
    return ind_dict

def update_index(ind_dict, spectrum, ind_val, redshift, y_scale):
    index_ax = ind_dict['ax']
    ind_dict_ = utils.shift_index(ind_dict, redshift)
    ind_dict_['text'].set_text(
        ind_dict_['name']+'='+str(np.round(ind_val,2)))
    if 'line' in ind_dict_.keys():
        continuum_flat = [j for sub in ind_dict_['continuum'] for j in sub]
        index_range = [min(continuum_flat), max(continuum_flat)]
        spec_range = np.where((spectrum[:,0] > index_range[0]) & (spectrum[:,0] < index_range[1]))
        ind_dict_['line'].set_ydata(spectrum[spec_range][:,1]*10**-y_scale)
        ind_dict_['line'].set_xdata(spectrum[spec_range][:,0])
        
        min_y = min(spectrum[spec_range][:,1]*10**-y_scale)
        max_y = max(spectrum[spec_range][:,1]*10**-y_scale)
        y_span = max_y - min_y
        index_ax.set_ylim([min_y-0.1*y_span, max_y+0.1*y_span])
        index_ax.set_xlim(index_range)
        
    if 'con_polys' in ind_dict_.keys():
        for i,poly in enumerate(ind_dict_['con_polys']):
            poly.set_verts([[[ind_dict_['continuum'][i][0],0],
                             [ind_dict_['continuum'][i][0],20],
                             [ind_dict_['continuum'][i][1],20],
                             [ind_dict_['continuum'][i][1],0]]])
            
    if 'feature_poly' in ind_dict_.keys():
        ind_dict_['feature_poly'].set_verts([[[ind_dict_['feature'][0],0],
                                              [ind_dict_['feature'][0],20],
                                              [ind_dict_['feature'][1],20],
                                              [ind_dict_['feature'][1],0]]])
        