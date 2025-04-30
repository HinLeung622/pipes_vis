import numpy as np

import bagpipes as pipes
from . import override_config
override_config.override_config(pipes)

from . import utils

# pe = plot element

def get_y_scale(spectrum, ymax=None):
    """ Sorts out y-axis scaling for plots with spectrum """
    if not ymax:
        ymax = 1.05*np.max(spectrum[:, 1])

    y_scale = int(np.log10(ymax))-1
    
    return y_scale,ymax
    

def add_bp_sfh(sfh, ax, color="black"):
    """ Adds a SFH black line in the given axis, in Bagpipes style """
    zvals=[0, 0.5, 1, 2, 4, 10]
    pe_sfh_line, = ax.plot(sfh[0], sfh[1],
            color=color, zorder=4, lw=2, alpha=1, ls="-", label=None)

    # Set limits
    ax.set_xlim([utils.cosmo.age(0).value, 0.])
    ax.set_ylim(bottom=0.)

    # Add redshift axis along the top
    pipes.plotting.add_z_axis(ax, zvals=zvals)

    # Add labels
    ax.set_ylabel("$\\mathrm{SFR\\ /\\ M_\\odot\\ \\mathrm{yr}^{-1}}$")
    ax.set_xlabel("$\\mathrm{Age\\ of\\ Universe\\ /\\ Gyr}$")
    return pe_sfh_line

def add_sfh_plot(params, sfh, input_logM, ax, sfh_color="black", sfh_label=None, 
                 z_line_color='red', z_line_label='obs redshift'):
    """
    Adds the elements in the sfh panel to the given axis,
    returns plot elements that can be later updated by sliders
    """
    z = params["redshift"]
    pe_sfh_line = add_bp_sfh(sfh, ax, color=sfh_color)
    age_at_z = utils.cosmo.age(z).value
    pe_z_line, = ax.plot([age_at_z,age_at_z], [0,ax.get_ylim()[1]*1000], color=z_line_color, 
                          label=z_line_label)
    pe_z_text = ax.annotate(str(np.round(z,3)), [age_at_z, 0.92*ax.get_ylim()[1]], color=z_line_color)
    ax.legend(loc='upper right')
    pe_input_logM_text = ax.text(0.985,0.88,r'pre-obs $\log_{10}M_*$ = '+str(np.round(input_logM,2)),
                   transform=ax.transAxes,fontsize=14, horizontalalignment='right',
                   verticalalignment='top', zorder=10)
    
    pe_bad_sfh_text = ax.text(0.5, 0.5, 'Oops! Bad SFH inputs, no SFH displayed', color='red',
                           fontsize=24, alpha=0, horizontalalignment='center', 
                           verticalalignment='center', transform=ax.transAxes, zorder=100)
    
    return pe_sfh_line, pe_z_line, pe_z_text, pe_input_logM_text, pe_bad_sfh_text

def add_bp_spectrum(spectrum, ax, zorder=4, z_non_zero=True, color="sandybrown",
                    ymax=None, lw=2., label=None, alpha=1, sub=False):
    """ 
    Add a spectrum to the passed axes. Adds errors if they are
    included in the spectrum object as a third column. 
    """
    y_scale, ymax = get_y_scale(spectrum, ymax)

    ax.set_ylim(0., ymax*10**-y_scale)
    ax.set_xlim(spectrum[0, 0], spectrum[-1, 0])

    # Plot the data
    pe_spec_line, = ax.plot(spectrum[:, 0], spectrum[:, 1]*10**-y_scale,
            color=color, zorder=zorder, lw=lw, label=label, alpha=alpha)

    # Sort out x tick locations
    pipes.plotting.auto_x_ticks(ax)
    if sub == True:
        ax.set_xscale('log')
        ax.set_yticks([])

    # Sort out axis labels.
    if sub == False:
        pipes.plotting.auto_axis_label(ax, y_scale, z_non_zero=z_non_zero)

    return y_scale, pe_spec_line

def add_main_spec(spectrum, ax, spec_lim, median_width=150, color="sandybrown", 
                  continuum_color='black', continuum_label='continuum'):
    """
    Adds the elements in the main spectrum plot to the given axis,
    returns plot elements that can be later updated by sliders and y scale of spectrum
    """
    run_med = utils.running_median(spectrum[:,0], spectrum[:,1], width=median_width)
    zoom_in_spec = spectrum[np.where((spectrum[:,0] >= spec_lim[0]) & 
                                     (spectrum[:,0] <= spec_lim[1]))]
    y_scale_spec,pe_spec_line = add_bp_spectrum(zoom_in_spec, ax, color=color)
    pe_run_med_line, = ax.plot(spectrum[:, 0], run_med*10**-y_scale_spec, color='k', alpha=0.3, zorder=5,
                            label='continuum')
    ax.legend(loc='upper right')
    pe_overflow_text = ax.text(0.5, 0.5, 'Fluxes too low, y-axis scale inaccurate!!', color='red',
                            fontsize=24, alpha=0, horizontalalignment='center', 
                            verticalalignment='center', transform=ax.transAxes, zorder=100)
    
    return pe_spec_line, pe_run_med_line, pe_overflow_text, y_scale_spec

def update_spec(spectrum, ax, pe_spec_line, sub=False, pe_overflow_text=None, change_xlims=False):
    """ updates a given spectrum line plot from a given ax with a given new spectrum """
    if pe_overflow_text is not None:
        pe_overflow_text.set_alpha(0.0)
    spec_ymax = 1.05*np.max(spectrum[:, 1])
    try:
        spec_y_scale = int(np.log10(spec_ymax))-1
        pe_spec_line.set_ydata(spectrum[:, 1]*10**-spec_y_scale)
        if change_xlims == True:
            pe_spec_line.set_xdata(spectrum[:, 0])
            ax.set_xlim(spectrum[0, 0], spectrum[-1, 0])
            pipes.plotting.auto_x_ticks(ax)
        ax.set_ylim(0., spec_ymax*10**-spec_y_scale)
        if sub == False:
            pipes.plotting.auto_axis_label(ax, spec_y_scale, z_non_zero=True)
        return spec_y_scale
    except OverflowError:
        if pe_overflow_text is not None:
            pe_overflow_text.set_alpha(1.0)

def add_residual(spectrum, ax, spec_lim, median_width=150, color="sandybrown"):
    """
    Adds the elements in the residual plot to the given axis,
    returns plot elements that can be later updated by sliders and y scale of spectrum
    """
    run_med = utils.running_median(spectrum[:,0], spectrum[:,1], width=median_width)
    residual = spectrum[:,1] / run_med
    pe_res_line, = ax.plot(spectrum[:, 0], residual, color=color, lw=2, zorder=1)
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
    in_range_res = residual[np.where((spectrum[:,0] >= spec_lim[0]) & 
                                     (spectrum[:,0] <= spec_lim[1]))]
    res_span = max(in_range_res) - min(in_range_res)
    ax.set_ylim([min(in_range_res)-0.1*res_span, max(in_range_res)+0.1*res_span])
    
    return pe_res_line

def update_residual(new_residual, pe_res_line, ax, model, spec_lim):
    pe_res_line.set_ydata(new_residual)
    in_range_res = new_residual[np.where((model.spectrum[:,0] >= spec_lim[0]) & 
                                         (model.spectrum[:,0] <= spec_lim[1]))]
    res_span = max(in_range_res) - min(in_range_res)
    ax.set_ylim([min(in_range_res)-0.1*res_span, max(in_range_res)+0.1*res_span])

def add_index_spectrum(ind_dict, spectrum, ind_val, redshift, y_scale=None, 
                       color_continuum='lightgray', color_feature='sandybrown', alpha=0.2):
    """ 
    Adds the elements in line index plots to a given axis, 
    returns plot elements that can later be updated by sliders and y scale of spectrum
    """
    if not y_scale:
        y_scale, ymax = get_y_scale(spectrum)
        
    # sort out units
    if ind_dict['type'] == 'break':
        units = ''
    else:
        if 'units' in ind_dict:
            if ind_dict['units'] == 'AA':
                units = '$\\AA$'
            elif ind_dict['units'] == 'mag':
                units = 'mag'
            else:
                units = ''
        else:
            units = '$\\AA$'
    
    index_ax = ind_dict['ax']
    ind_dict_ = utils.shift_index(ind_dict, redshift)
    if ind_dict_['type'] in ['EW', 'break']:
        continuum_flat = [j for sub in ind_dict_['continuum'] for j in sub]
        index_range = [min(continuum_flat), max(continuum_flat)]
        spec_range = np.where((spectrum[:,0] > index_range[0]) & (spectrum[:,0] < index_range[1]))
        pe_ind_line = index_ax.plot(spectrum[spec_range][:,0], spectrum[spec_range][:,1]*10**-y_scale, color='k')
        ind_dict['line'] = pe_ind_line[0]
        
        pe_ind_y_scale_text = index_ax.text(0.0,1.01, r'$\times 10^{'+str(y_scale)+r'}$', 
                                         ha='left', va='bottom', transform=index_ax.transAxes)
        ind_dict['y_scale_text'] = pe_ind_y_scale_text

        ylims = index_ax.get_ylim()
        pe_con_polys = []
        for con in ind_dict_['continuum']:
            pe_con_poly = index_ax.fill_between(con, [0,0], [20,20], color='lightgray', alpha=0.2)
            pe_con_polys.append(pe_con_poly)
        ind_dict['con_polys'] = pe_con_polys
        if ind_dict_['type'] == 'EW':
            pe_feature_poly = index_ax.fill_between(ind_dict_['feature'], [0,0], [20,20], color='sandybrown', alpha=0.2)
            ind_dict['feature_poly'] = pe_feature_poly
        index_ax.set_ylim(ylims)
        index_ax.set_xlim(index_range)
        index_ax.tick_params(direction="in")
        pe_ind_text = index_ax.text(1.0, 1.01, ind_dict_['name']+'='+str(np.round(ind_val,2))+units, 
                                 ha='right', va='bottom', transform=index_ax.transAxes)
    else:
        index_ax.set_xticks([])
        pe_ind_text = index_ax.text(0.5, 0.5, ind_dict_['name']+'='+str(np.round(ind_val,2))+units, 
                                 ha='center', va='center', transform=index_ax.transAxes)
    
    ind_dict['text'] = pe_ind_text
    
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
        
    if 'y_scale_text' in ind_dict_.keys():
        ind_dict_['y_scale_text'].set_text(r'$\times 10^{'+str(y_scale)+r'}$')
        
    if len(index_ax.get_ylabel()) > 0:
        index_ax.yaxis.set_label("$\\mathrm{f_{\\lambda}}\\ \\mathrm{/\\ 10^{"
                            + str(y_scale)
                            + "}\\ erg\\ s^{-1}\\ cm^{-2}\\ \\AA^{-1}}$")
        
def update_zmet_plot(ax, pe_zmet_line, zmet_evo):
    pe_zmet_line.set_ydata(zmet_evo)
    max_y = np.nanmax(zmet_evo)
    ax.set_ylim(top = 1.1*max_y)

def add_input_photometry(galaxy, ax, zorder=6, y_scale=None, ptsize=40, lw=1.,
                         marker="o", color="blue"):
    """ Adds photometric data to the passed axes without doing any
    manipulation of the axes or labels. 
    (Adopted from bagpipes.plotting.add_input_photometry_linear)"""

    photometry = np.copy(galaxy.photometry)

    # Plot the data
    pe_input_phot_errbar = ax.errorbar(photometry[:, 0], photometry[:, 1]*10**-y_scale,
        yerr=photometry[:, 2]*10**-y_scale, lw=lw,
        linestyle=" ", capsize=3, capthick=lw, zorder=zorder-1,
        color="black")

    pe_input_phot = ax.scatter(photometry[:, 0], photometry[:, 1]*10**-y_scale, color=color,
        s=ptsize, zorder=zorder, linewidth=lw, facecolor=color,
        edgecolor="black", marker=marker)

    return pe_input_phot, pe_input_phot_errbar

def update_input_photometry(galaxy, pe_input_phot, pe_input_phot_errbar, y_scale):
    photometry = np.copy(galaxy.photometry)

    new_offset = np.vstack([photometry[:, 0], photometry[:, 1]*10**-y_scale]).T
    pe_input_phot.set_offsets(new_offset)

    # structure of errorbar object is line, (bottom_caps, top_caps), vertical_bars
    new_yerr = photometry[:, 2]*10**-y_scale
    pe_input_phot_errbar[0].set_ydata(new_offset[:,1])
    pe_input_phot_errbar[1][0].set_ydata(new_offset[:,1] - new_yerr)
    pe_input_phot_errbar[1][1].set_ydata(new_offset[:,1] + new_yerr)
    np.vstack([photometry[:, 0], new_offset[:,1] - new_yerr])
    # structure of segments is [[[x0, y0_lower], [x0, y0_higher]], [[x1, y1_lower], [x1, y1_higher]]...]
    new_segments = np.stack(
        [np.vstack([photometry[:, 0], new_offset[:,1] - new_yerr]).T, 
         np.vstack([photometry[:, 0], new_offset[:,1] + new_yerr]).T]
        ).transpose(1,0,2)
    pe_input_phot_errbar[2][0].set_segments(new_segments)

def add_model_photometry(eff_wavs, photometry, ax, y_scale=None, color="sandybrown", s=100, zorder=5):
    """ Adds model photometry to the passed axis. """

    pe_model_phot = ax.scatter(
        eff_wavs, photometry*10**-y_scale, facecolor=color, edgecolor='k', s=s, zorder=zorder, marker='s'
        )
    
    return pe_model_phot

def update_model_photometry(eff_wavs, new_photometry, pe_model_phot, y_scale):
    photometry_ = np.copy(new_photometry)
    new_offset = np.vstack([eff_wavs, photometry_*10**-y_scale]).T
    pe_model_phot.set_offsets(new_offset)
