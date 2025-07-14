import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox, CheckButtons
import copy
import os
#from time import perfcounter

import bagpipes as pipes
from . import override_config
override_config.override_config(pipes)

from . import utils
from . import plotting
from . import slider_params


class visualizer:
    """
    Object that holds and does the under-the-hood work for showing the galaxy SFH/parameters
    -> spectrum GUI visualizer.
    """
    def __init__(self, init_components, wavelengths=None,
                 load_data_func=None, spec_units='ergscma', phot_units='mujy', 
                 spectrum_exists=True, photometry_exists=True, 
                 filt_list=None, 
                 index_list=None):
        self.init_comp = init_components
        self.input_galaxy_kwarg = {
            "load_data": load_data_func,
            "spec_units": spec_units,
            "phot_units": phot_units,
            "spectrum_exists": spectrum_exists,
            "photometry_exists": photometry_exists,
            "filt_list": filt_list
        }
        self.init_spec_lim = self.init_comp['spec_lim'].copy()
        self.index_list = index_list
        self.GUI_initialized = False
        self.additional_plots = False
        # full spectrum wavelengths
        if wavelengths is None:
            wavelengths_fine = np.linspace(1000,10000, 10000)
            wavelengths_coarse = np.linspace(10000,50000,2000)[1:]
            self.wavelengths = np.concatenate([wavelengths_fine, wavelengths_coarse])
        else:
            self.wavelengths = wavelengths
        # colours used for the various plot lines
        self.plot_colors = {'sfh':"black", 'z':"red", 'spectrum':"sandybrown", 
                            'continuum':"black", 'zoom':"black", 'photometry':"red",
                            'input_photometry':'blue', 'input_spec':'dodgerblue',
                            "index_continuum":"lightgray", "index_feature":"sandybrown"}
        # width of spectrum (in AA) considered when calculating running median
        self.median_width = 150
        # adjustment values that determine sizes and arrangement of sliders
        self.left_x = 0.12
        self.right_x = 0.57
        self.max_y = 0.28      # y of top most slider
        self.step_y = 0.03     # y gap between consecutive sliders
        self.width = 0.33
        self.height = 0.025
        self.slider_gap = 0.07 # x gap between columns of sliders
        self.max_n_rows = 9
        self.textbox_gap = 0.03
        # other adjustment values for other elements of the plot
        self.sub_ax_arg = [0.7,0.05,0.28,0.4]          # location and scaling of full spectrum subplot
        self.bottom_adjust_val = 0.36                  # scaling value used on all plots to make space
        self.top_adjust_val = 0.95                     # scaling value used on all plots to make space
                                                       # for the sliders
        self.slider_colors = {'sfh':['cornflowerblue','orange','lime','tomato','darkorchid'],
                              'misc':'gold', 'reset_button':'lightgoldenrodyellow'}       # colour of the sliders
        self.reset_button_arg = [0.8, 0.025, 0.1, 0.04]           # location and scaling of the reset button
        self.moreplots_check_arg = [0.7, 0.025, 0.04,0.04]        # location and scaling of checkbox for more plots
        
        self.plot_input_phot = False
        self.plot_input_spec = False
        if load_data_func is not None:
            self.input_galaxy = pipes.galaxy('_', **self.input_galaxy_kwarg)
            self.plot_input_phot = photometry_exists
            self.plot_input_spec = spectrum_exists

    def static_plot(self, show=True, figsize=(13,9)):
        """ 
        Creates the figure, lines, texts and annotations. Returns figure and axes (in a list)
        for further manipulation by the user if needed.
        """
        self.fig = plt.figure(figsize=figsize)

        if self.index_list is None:
            gs1 = matplotlib.gridspec.GridSpec(13, 1, hspace=0., wspace=0.)
        else:
            index_ncols = -(-len(self.index_list)//5)
            gs0 = self.fig.add_gridspec(1, 3+index_ncols)
            gs1 = gs0[:3].subgridspec(13, 1, hspace=0., wspace=0.)
            gs2 = gs0[3:].subgridspec(5, index_ncols, hspace=0.5, wspace=0.2)
            # indices plots
            for i in range(len(self.index_list)):
                self.index_list[i]['ax'] = plt.subplot(gs2[i%5,i//5])
            
        self.ax1 = plt.subplot(gs1[0:5])        #SFH plot
        self.ax2 = plt.subplot(gs1[6:11])       #main spectrum plot
        self.ax3 = plt.subplot(gs1[11:])        #residual plot
                

        init_input_logM, self.total_sfh, init_custom_sfh = utils.create_sfh(self.init_comp)

        self.pe_sfh_line, self.pe_z_line, self.pe_z_text, self.pe_input_logM_text, self.pe_bad_sfh_text \
            = plotting.add_sfh_plot(self.init_comp, self.total_sfh, init_input_logM, self.ax1,
                                    sfh_color=self.plot_colors['sfh'], z_line_color=self.plot_colors['z'])

        self.model = pipes.model_galaxy(utils.make_pipes_components(self.init_comp, init_custom_sfh),
                                        spec_wavs=self.wavelengths, filt_list=self.input_galaxy_kwarg['filt_list'])

        # full spectrum in inset
        self.sub_ax = plt.axes([0,0,1,1])
        self._adjust_subax_pos()

        sub_y_scale_spec,self.pe_sub_spec_line = plotting.add_bp_spectrum(self.model.spectrum, self.sub_ax, sub=True,
                                                                       color=self.plot_colors['spectrum'])
        self.pe_spec_zoom_poly = self.sub_ax.fill_between(self.init_spec_lim, [0]*2, [20]*2, color=self.plot_colors['zoom'], 
                                                          alpha=0.1)

        # the main spectrum plot
        self.pe_spec_line, self.pe_run_med_line, self.pe_overflow_text, y_scale_spec \
            = plotting.add_main_spec(self.model.spectrum, self.ax2, self.init_spec_lim, median_width=self.median_width, 
                                     color=self.plot_colors['spectrum'],
                                     continuum_color=self.plot_colors['continuum'])

        # the residual plot
        self.pe_res_line = plotting.add_residual(self.model.spectrum, self.ax3, self.init_spec_lim, median_width=self.median_width, 
                                                 color=self.pe_spec_line.get_color())
        
        # indices plots
        if self.index_list is not None:
            self.index_names = [ind["name"] for ind in self.index_list]

            # calculate and plot indices
            self.indices = np.zeros(len(self.index_list))
            for i in range(self.indices.shape[0]):
                self.indices[i] = pipes.input.spectral_indices.measure_index(
                    self.index_list[i], self.model.spectrum, self.init_comp["redshift"])
                
                self.index_list[i] = plotting.add_index_spectrum(
                    self.index_list[i], self.model.spectrum, self.indices[i], self.init_comp["redshift"], 
                    y_scale=y_scale_spec, color_continuum=self.plot_colors['index_continuum'],
                    color_feature=self.plot_colors['index_feature'], alpha=0.2
                    )
                
        if self.plot_input_phot:
            # add input photometry and spectrum
            self.pe_input_phot_scatter, self.pe_input_phot_errbar = plotting.add_input_photometry(
                self.input_galaxy, self.ax2, y_scale=y_scale_spec, 
                color=self.plot_colors['input_photometry']
                )
            self.pe_sub_input_phot_scatter, self.pe_sub_input_phot_errbar = plotting.add_input_photometry(
                self.input_galaxy, self.sub_ax, y_scale=sub_y_scale_spec, 
                color=self.plot_colors['input_photometry'], ptsize=10
                )
            self.sub_ax.set_ylabel(None)

            # add predicted photometry
            self.pe_model_phot = plotting.add_model_photometry(
                self.model.filter_set.eff_wavs, self.model.photometry, self.ax2, y_scale=y_scale_spec, color=self.plot_colors['photometry']
                )
            self.pe_sub_model_phot = plotting.add_model_photometry(
                self.model.filter_set.eff_wavs, self.model.photometry, self.sub_ax, y_scale=sub_y_scale_spec, color=self.plot_colors['photometry'], s=25
                )

        #plt.tight_layout()
        
        if show:
            plt.show()
            
        return self.fig, [self.ax1, self.ax2, self.ax3, self.sub_ax]
    
    def index_static_plot(self, show=True, figsize=(13,9)):
        """ 
        Indices only
        Creates the figure, lines, texts and annotations. Returns figure and axes (in a list)
        for further manipulation by the user if needed.
        """
        self.fig = plt.figure(figsize=figsize)

        gs1 = matplotlib.gridspec.GridSpec(13, 1, hspace=0., wspace=0.)

        self.ax1 = plt.subplot(gs1[0:5])        #SFH plot
        # indices plots
        index_ncols = -(-len(self.index_list)//2)
        gs2 = gs1[6:].subgridspec(2, index_ncols, hspace=0.5, wspace=0.2)
        for i in range(len(self.index_list)):
            self.index_list[i]['ax'] = plt.subplot(gs2[i%2,i//2])


        init_input_logM, self.total_sfh, init_custom_sfh = utils.create_sfh(self.init_comp)

        self.pe_sfh_line, self.pe_z_line, self.pe_z_text, self.pe_input_logM_text, self.pe_bad_sfh_text \
            = plotting.add_sfh_plot(self.init_comp, self.total_sfh, init_input_logM, self.ax1,
                                    sfh_color=self.plot_colors['sfh'], z_line_color=self.plot_colors['z'])

        self.model = pipes.model_galaxy(utils.make_pipes_components(self.init_comp, init_custom_sfh),
                                        spec_wavs=self.wavelengths)

        # indices plots
        self.index_names = [ind["name"] for ind in self.index_list]

        # calculate and plot indices
        y_scale_spec,ymax = plotting.get_y_scale(self.model.spectrum)
        self.indices = np.zeros(len(self.index_list))
        for i in range(self.indices.shape[0]):
            self.indices[i] = pipes.input.spectral_indices.measure_index(
                self.index_list[i], self.model.spectrum, self.init_comp["redshift"])

            self.index_list[i] = plotting.add_index_spectrum(
                self.index_list[i], self.model.spectrum, self.indices[i], self.init_comp["redshift"], 
                y_scale=y_scale_spec,
                color_continuum=self.plot_colors['index_continuum'],
                color_feature=self.plot_colors['index_feature'], alpha=0.2
                )

        # add labels to index plots
        pipes.plotting.auto_axis_label(self.index_list[0]['ax'], y_scale_spec)
        self.index_list[0]['ax'].yaxis.set_label_coords(-0.05, -0.2)

        #plt.tight_layout()

        if show:
            plt.show()

        return self.fig
        
    def GUI(self, figsize=(13,12), index_only=False):
        """ 
        Creates the figure, lines, texts and annotations that will be manipulated, and also the 
        interactive elements.
        """
        # check if the GUI has already been initialized. If it is, just show the figure from the previous state
        if self.GUI_initialized:
            self.fig.show()
        else:
            # initializing GUI plot
            if index_only:
                self.index_static_plot(show=False, figsize=figsize)
            else:
                self.static_plot(show=False, figsize=figsize)
            self.spec_lim = self.init_spec_lim.copy()
    
            # adjust the main plots to make room for the sliders
            plt.subplots_adjust(bottom=self.bottom_adjust_val, top=self.top_adjust_val) 
            self._adjust_subax_pos()
    
            # make sliders
            self.sliders_ax, self.sliders, current_right_column_pos = self.make_sliders()
    
            # make input textboxes for wavelength limits
            if current_right_column_pos[0] == self.right_x:
                right_edge = self.reset_button_arg[0] - 0.03
            else:
                right_edge = self.right_x+self.width
            textbox_width = (right_edge-current_right_column_pos[0]-self.textbox_gap)/2
    
            self.ax_spec_min = plt.axes([current_right_column_pos[0], 
                                         current_right_column_pos[1], textbox_width, self.height])
            self.spec_min_box = TextBox(self.ax_spec_min, r'$\lambda_{min}$', initial=self.init_spec_lim[0])
            self.spec_min_box.on_submit(self.submit_min)
    
            self.ax_spec_max = plt.axes([current_right_column_pos[0]+textbox_width+self.textbox_gap, 
                                         current_right_column_pos[1], textbox_width, self.height])
            self.spec_max_box = TextBox(self.ax_spec_max, r'$\lambda_{max}$', initial=self.init_spec_lim[1])
            self.spec_max_box.on_submit(self.submit_max)
    
            # register the update function with each slider
            for key in self.sliders.keys():
                self.sliders[key].on_changed(self.update)
    
            # Create a `matplotlib.widgets.Button` to reset all sliders to initial values.
            self.resetax = plt.axes(self.reset_button_arg)
            self.reset_button = Button(self.resetax, 'Reset', color=self.slider_colors['reset_button'], hovercolor='0.975')
    
            self.reset_button.on_clicked(self.reset)
    
            # Create a `matplotlib.widgets.CheckButton` to toggle show or not additional plots
            self.moreplotsax = plt.axes(self.moreplots_check_arg)
            self.moreplots_check = CheckButtons(self.moreplotsax, ['additional plots'])
    
            self.moreplots_check.on_clicked(self.toggle_additional_plots)
    
            # additional plots toggle bool
            self.reset_additional_plots()
    
            plt.show()
    
            self.GUI_initialized = True
        
    def make_one_slider(self, x_pos, y_pos, width, height, label, lims, init_val, color):
        """ makes a single slider """
        bg_color = list(matplotlib.colors.to_rgba(color))
        bg_color[3] = 0.25
        ax_slider = plt.axes([x_pos, y_pos, width, height], facecolor=bg_color, alpha=0.5)

        if int(matplotlib.__version__.split('.')[1]) >= 5:
            # new slider format
            the_slider = Slider(
                ax=ax_slider,
                label=label,
                valmin=lims[0],
                valmax=lims[1],
                valinit=init_val,
                color=color,
                track_color=bg_color
            )
        else:
            # old slider format
            the_slider = Slider(
                ax=ax_slider,
                label=label,
                valmin=lims[0],
                valmax=lims[1],
                valinit=init_val,
                color=color
            )

        return ax_slider, the_slider

    def make_sliders(self):
        """
        detects the lists of sliders that has to be made from the initial components, 
        then uses the slider_params.slider_lib to arrange and create the proper list of sliders at the correct locations,
        with min/max values either from the lib or user-specified values.
        Returns the sliders dictionary and the "upper-most empty y coord" on both slider columns
        """
        # loop through the components to 1. identify which components are reqiured to create sliders for
        # 2. compile info such as initial values and custom slider limits into a list of the same order
        required_components = []
        init_vals = []
        custom_extremes = {}
        slider_names = []
        repeat_info = {}    # element 0 holds the name of the repeat, 1 holds the number of repeats
        sfh_types = dir(pipes.models.star_formation_history)
        for key in self.init_comp.keys():
            # build up repeat_info
            if key in sfh_types or key[:-1] in sfh_types:
                if key[:-1] in sfh_types:
                    if key[:-1] not in repeat_info.keys():
                        repeat_info[key[:-1]] = 1
                    else:
                        repeat_info[key[:-1]] += 1
                        
                sfh_dict = self.init_comp[key]
                for sfh_key in sfh_dict.keys():
                    if sfh_key[-5:] != '_lims':
                        if key in sfh_types and key+':'+sfh_key in slider_params.slider_lib.keys():
                            required_components.append(key+':'+sfh_key)
                            slider_names.append(key+':'+sfh_key)
                            init_vals.append(sfh_dict[sfh_key])
                        elif key[:-1] in sfh_types and key[:-1]+':'+sfh_key in slider_params.slider_lib.keys():
                            required_components.append(key[:-1]+':'+sfh_key)
                            slider_names.append(key+':'+sfh_key)
                            init_vals.append(sfh_dict[sfh_key])
                        if sfh_key+'_lims' in sfh_dict.keys():
                            custom_extremes[key+':'+sfh_key] = sfh_dict[sfh_key+'_lims']
                            
            elif key in ['dust', 'nebular']:
                for sub_key in self.init_comp[key].keys():
                    if sub_key[-5:] != '_lims':
                        if sub_key != 'type' and key+':'+sub_key in slider_params.slider_lib.keys():
                            required_components.append(key+':'+sub_key)
                            slider_names.append(key+':'+sub_key)
                            init_vals.append(self.init_comp[key][sub_key])
                            if sub_key+'_lims' in self.init_comp[key].keys():
                                custom_extremes[key+':'+sub_key] = self.init_comp[key][sub_key+'_lims']
                                
            elif key != 'spec_lim' and key[-5:] != '_lims' and key in slider_params.slider_lib.keys():
                required_components.append(key)
                slider_names.append(key)
                init_vals.append(self.init_comp[key])
                if key+'_lims' in self.init_comp.keys():
                    custom_extremes[key] = self.init_comp[key+'_lims']
            
        # loop through the reqiured components a second time to 1. seperate left and right sliders
        # 2. take according standard slider properties from slider_params.slider_lib
        # 3. replace specified slider limits with custom limits
        # 4. sort out ordering with repeated sfh components
        left_components = []
        left_priorities = []
        right_components = []
        right_priorities = []
        for i,component in enumerate(required_components):
            c_dict = slider_params.slider_lib[component].copy()
            c_dict['name'] = slider_names[i]
            c_dict['init_val'] = init_vals[i]
            
            if ':' in component:
                if component.split(':')[0] in repeat_info.keys():
                    # found a parameter of a repeated sfh component
                    sfh_type = slider_names[i].split(':')[0][:-1]
                    rep_order = int(slider_names[i].split(':')[0][-1])
                    c_dict['priority'] = slider_params.sfh_priorities[sfh_type] + \
                        (rep_order-1)*10/repeat_info[sfh_type] + \
                        (c_dict['priority']-slider_params.sfh_priorities[sfh_type])/repeat_info[sfh_type]
            
            if slider_names[i] in custom_extremes.keys():
                c_dict['lims'] = custom_extremes[slider_names[i]]
            if c_dict['side'] == 'left':
                left_priorities.append(c_dict['priority'])
                left_components.append(c_dict)
            elif c_dict['side'] == 'right':
                right_priorities.append(c_dict['priority'])
                right_components.append(c_dict)

        left_components = [left_components[i] for i in np.argsort(left_priorities)]
        right_components = [right_components[i] for i in np.argsort(right_priorities)]

        # now actually making the sliders with the info ordered above
        sliders_ax = {}
        sliders = {}
        sfh_comp_list = []
        
        # left columns, includes auto ordering in cases with more sliders than self.max_n_rows,
        # and auto assigning of colours to groups in SFH components
        left_n_cols = np.ceil(len(left_components)/self.max_n_rows)
        left_width = (self.width-(left_n_cols-1)*self.slider_gap)/left_n_cols
        left_row_no = np.arange(len(left_components))%self.max_n_rows
        left_col_no = np.floor(np.arange(len(left_components))/self.max_n_rows)
        for i,c_dict in enumerate(left_components):
            slider_name = c_dict['name']
            # check for each group of SFH params and give them sorted slider background colours accordingly
            if ':' in slider_name:
                if slider_name.split(':')[0] not in sfh_comp_list:
                    sfh_comp_list.append(slider_name.split(':')[0])
                bg_color = self.slider_colors['sfh'][(len(sfh_comp_list)-1)%len(self.slider_colors['sfh'])]
            else:
                bg_color = self.slider_colors['misc']
            x_pos = left_col_no[i]*(left_width+self.slider_gap) + self.left_x
            y_pos = self.max_y - left_row_no[i]*self.step_y
            sliders_ax[slider_name], sliders[slider_name] \
                = self.make_one_slider(x_pos, y_pos, left_width, self.height, c_dict['label'], c_dict['lims'],
                                    c_dict['init_val'], bg_color)
        
        # right columns, includes auto ordering in cases with more sliders than self.max_n_rows,
        # but always the same colour
        right_n_cols = np.ceil((len(right_components)+1)/self.max_n_rows)
        right_width = (self.width-(right_n_cols-1)*self.slider_gap)/right_n_cols
        right_row_no = np.arange(len(right_components)+1)%self.max_n_rows
        right_col_no = np.floor(np.arange(len(right_components)+1)/self.max_n_rows)
        for i,c_dict in enumerate(right_components):
            slider_name = c_dict['name']
            bg_color = self.slider_colors['misc']
            x_pos = right_col_no[i]*(right_width+self.slider_gap) + self.right_x
            y_pos = self.max_y - right_row_no[i]*self.step_y
            sliders_ax[slider_name], sliders[slider_name] \
                = self.make_one_slider(x_pos, y_pos, right_width, self.height, c_dict['label'], c_dict['lims'],
                                    c_dict['init_val'], bg_color)
        
        # calculate the coords of the next empty slot for putting the text boxes
        current_right_column_pos = [right_col_no[-1]*(right_width+self.slider_gap) + self.right_x,
                                    self.max_y - right_row_no[-1]*self.step_y]
            
        return sliders_ax, sliders, current_right_column_pos

    def update_component_dict(self):
        """ create updated components dictionary """
        sfh_dict_list = []
        sfh_types = dir(pipes.models.star_formation_history)
        for key in self.init_comp.keys():
            if key in sfh_types or key[:-1] in sfh_types:
                sfh_dict = {}
                for sfh_key in self.init_comp[key].keys():
                    if sfh_key[-5:] == '_lims':
                        # first throw away the limits
                        continue
                    elif key+':'+sfh_key in self.sliders.keys():
                        sfh_dict[sfh_key] = self.sliders[key+':'+sfh_key].val
                    else:
                        # a catch all for anything that is not modifiable through sliders
                        sfh_dict[sfh_key] = self.init_comp[key][sfh_key]
                sfh_dict_list.append(sfh_dict)
        
        other_modules = {}
        for module in ['dust', 'nebular']:
            if module in self.init_comp.keys():
                module_dict = {}
                for key in self.init_comp[module].keys():
                    if key[-5:] == '_lims':
                        # first throw away the limits
                        continue
                    elif key == 'type':
                        module_dict["type"] = self.init_comp[module]['type']
                    elif module+':'+key in self.sliders.keys():
                        module_dict[key] = self.sliders[module+':'+key].val
                    else:
                        # a catch all for anything that is not modifiable through sliders
                        module_dict[key] = self.init_comp[module][key]
                other_modules[module] = module_dict

        self.new_comp = {}
        sfh_index = 0
        for key in self.init_comp.keys():
            if key == 'spec_lim' and key[-5:] == '_lims':
                # first throw away the limits
                continue
            elif key in sfh_types or key[:-1] in sfh_types:
                self.new_comp[key] = sfh_dict_list[sfh_index]
                sfh_index += 1
            elif key == 'dust':
                self.new_comp["dust"] = other_modules['dust']
            elif key == 'nebular':
                self.new_comp["nebular"] = other_modules['nebular']
            elif key in self.sliders.keys():
                self.new_comp[key] = self.sliders[key].val
            else:
                # a catch all for anything that is not modifiable through sliders
                self.new_comp[key] = self.init_comp[key]

    def update(self, val):
        """ The function to be called anytime a slider's value changes """
        # get new_comp going
        self.update_component_dict()

        #update sfh plot
        age_at_z = utils.cosmo.age(self.new_comp["redshift"]).value
        self.pe_z_line.set_xdata([age_at_z,age_at_z])
        input_logM, self.total_sfh, custom_sfh = utils.create_sfh(self.new_comp)
        self.pe_sfh_line.set_ydata(self.total_sfh[1])
        if max(self.total_sfh[1])>0:
            self.ax1.set_ylim(top=1.05*max(self.total_sfh[1]))
            self.pe_bad_sfh_text.set_alpha(0.0)
        else:
            self.ax1.set_ylim(top=1)
            self.pe_bad_sfh_text.set_alpha(1.0)
        self.pe_z_text.remove()
        self.pe_z_text = self.ax1.annotate(str(np.round(self.new_comp["redshift"],3)), 
                                        [age_at_z, 0.92*self.ax1.get_ylim()[1]], color=self.pe_z_line.get_color())
        self.pe_input_logM_text.set_text('pre-obs log10M = '+str(np.round(input_logM,2)))

        #update model components
        self.model.update(utils.make_pipes_components(self.new_comp, custom_sfh))
        
        run_med = utils.running_median(self.model.spectrum[:,0], self.model.spectrum[:,1], 
                                       width=self.median_width)
        residual = self.model.spectrum[:,1] / run_med

        #update subplot full spectrum
        if hasattr(self, 'sub_ax'):
            sub_y_scale_spec = plotting.update_spec(self.model.spectrum, self.sub_ax, self.pe_sub_spec_line, sub=True)

            if self.plot_input_phot:
                # update input galaxy data
                plotting.update_input_photometry(
                    self.input_galaxy, self.pe_sub_input_phot_scatter, self.pe_sub_input_phot_errbar, sub_y_scale_spec
                    )
                
                # update model photometry
                plotting.update_model_photometry(
                    self.model.filter_set.eff_wavs, self.model.photometry, self.pe_sub_model_phot, sub_y_scale_spec
                )

        #update main spectrum plot
        if hasattr(self, 'ax2'):
            zoom_in_spec = self.model.spectrum[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                                        (self.model.spectrum[:,0] <= self.spec_lim[1]))]
            y_scale_spec = plotting.update_spec(zoom_in_spec, self.ax2, self.pe_spec_line, pe_overflow_text=self.pe_overflow_text)
            if y_scale_spec is not None:
                # update running median line
                self.pe_run_med_line.set_ydata(run_med*10**-y_scale_spec)

                #update residual plot
                plotting.update_residual(residual, self.pe_res_line, self.ax3, self.model, self.spec_lim)

                if self.plot_input_phot:
                    # update input galaxy data
                    plotting.update_input_photometry(
                        self.input_galaxy, self.pe_input_phot_scatter, self.pe_input_phot_errbar, y_scale_spec
                        )
                    
                    # update model photometry
                    plotting.update_model_photometry(
                        self.model.filter_set.eff_wavs, self.model.photometry, self.pe_model_phot, y_scale_spec
                    )
            
        #update indices
        if self.index_list is not None:
            if hasattr(self, 'ax2') == False:
                y_scale_spec,ymax = plotting.get_y_scale(self.model.spectrum)
            for i in range(len(self.index_list)):
                self.indices[i] = pipes.input.spectral_indices.measure_index(
                    self.index_list[i], self.model.spectrum, self.new_comp["redshift"])
                plotting.update_index(self.index_list[i], self.model.spectrum, self.indices[i], 
                                      self.new_comp["redshift"], y_scale_spec)
                        
        if self.additional_plots:
            # mass weighted metallicity plot
            self.zmet_evo = utils.get_ceh(self.total_sfh[0], self.new_comp, self.model)
            plotting.update_zmet_plot(self.ax_zmet, self.zmet_line, self.zmet_evo)
        
        self.fig.canvas.draw_idle()
        
    def get_model(self):
        """ Returns the bagpipes model galaxy at its current state """
        return self.model
    
    def _adjust_subax_pos(self):
        # calculate plot coordinates relative to ax2
        if hasattr(self, 'sub_ax') and hasattr(self, 'ax2'):
            l, b, w, h = self.ax2.get_position().bounds
            sub_ax_arg_transformed = [
                self.sub_ax_arg[0] * w + l,
                self.sub_ax_arg[1] * h + b,
                self.sub_ax_arg[2] * w,
                self.sub_ax_arg[3] * h,
            ]
            self.sub_ax.set_position(sub_ax_arg_transformed)

    def get_spectrum(self, full=False):
        """ Returns the spectrum arrray at its current state.
        If full==True, returns the whole wavelength range.
        Column 0 = wavelength in Anstroms
        Column 1 = flux in erg s^-1 cm^-2 AA^-1 """
        if full:
            return self.model.spectrum
        else:
            return self.model.spectrum[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                                (self.model.spectrum[:,0] <= self.spec_lim[1]))]

    def submit_min(self, text):
        """
        adjust the main spectrum plot's limits and the full spec's shaded region
        on submit of a new wavelength lower limit value
        """
        self.spec_lim[0] = eval(text)
        zoom_in_spec = self.model.spectrum[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                               (self.model.spectrum[:,0] <= self.spec_lim[1]))]
        plotting.update_spec(zoom_in_spec, self.ax2, self.pe_spec_line, change_xlims=True)
        self.ax3.set_xlim(self.spec_lim)
        pipes.plotting.auto_x_ticks(self.ax3)
        self.pe_spec_zoom_poly.set_verts([[[self.spec_lim[0],0],
                                           [self.spec_lim[0],20],
                                           [self.spec_lim[1],20],
                                           [self.spec_lim[1],0]]])
        self.fig.canvas.draw_idle()

    def submit_max(self, text):
        """
        adjust the main spectrum plot's limits and the full spec's shaded region
        on submit of a new wavelength upper limit value
        """
        self.spec_lim[1] = eval(text)
        zoom_in_spec = self.model.spectrum[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                               (self.model.spectrum[:,0] <= self.spec_lim[1]))]
        plotting.update_spec(zoom_in_spec, self.ax2, self.pe_spec_line, change_xlims=True)
        self.ax3.set_xlim(self.spec_lim)
        pipes.plotting.auto_x_ticks(self.ax3)
        self.pe_spec_zoom_poly.set_verts([[[self.spec_lim[0],0],
                                           [self.spec_lim[0],20],
                                           [self.spec_lim[1],20],
                                           [self.spec_lim[1],0]]])
        self.fig.canvas.draw_idle()
        
    def reset(self, event):
        """ resets everything to the default states """
        for key in self.sliders.keys():
            self.sliders[key].reset()
        self.spec_min_box.set_val(self.init_spec_lim[0])
        self.spec_max_box.set_val(self.init_spec_lim[1])
        
    def toggle_additional_plots(self, label):
        """
        Toggling the display of additional plots by the check box. Currently
        plots include:
            Mass weighted metallicity plot, competible with time-varying models,
            plotted with SFH plot panel
        """
        if self.additional_plots:
            self.additional_plots = False     #toggle
            self.ax_zmet.set_axis_off()
            self.zmet_line.set_visible(False)
        else:
            self.additional_plots = True     #toggle
            if hasattr(self, 'ax_zmet'):
                self.zmet_evo = utils.get_ceh(self.total_sfh[0], self.new_comp, self.model)
                plotting.update_zmet_plot(self.ax_zmet, self.zmet_line, self.zmet_evo)
                self.ax_zmet.set_axis_on()
                self.zmet_line.set_visible(True)
            else:
                if hasattr(self, 'new_comp'):
                    self.zmet_evo = utils.get_ceh(self.total_sfh[0], self.new_comp, self.model)
                else:
                    self.zmet_evo = utils.get_ceh(self.total_sfh[0], self.init_comp, self.model)
                self.ax_zmet = self.ax1.twinx()
                line = self.ax_zmet.plot(self.total_sfh[0], self.zmet_evo, color='sandybrown')
                self.zmet_line = line[0]
                self.ax_zmet.set_ylim([0,1.1*np.nanmax(self.zmet_evo)])
                self.ax_zmet.set_ylabel(r'$Z_*/Z_\odot$', color='sandybrown')
                self.ax_zmet.tick_params(axis='y', labelcolor='sandybrown')
                self.ax_zmet.tick_params(direction="in")
            
        self.fig.canvas.draw_idle()
    
    def reset_additional_plots(self):
        """ Resets the additional plots related variables """
        self.additional_plots = False
        if hasattr(self, 'ax_zmet'):
            delattr(self, 'ax_zmet')
            
    def quick_lick(self, index_plot=True):
        """
        Quick access to a GUI version that shows a list of the lick indices
        """
        # set up index_list
        vis_path = os.path.dirname(os.path.realpath(__file__))
        lick_index = np.loadtxt(vis_path+'/lick_index.txt',usecols=[1,2,3,4,5,6],comments='#')
        names = np.loadtxt(vis_path+'/lick_index.txt',usecols=[0],dtype=str,comments='#')
        units = np.loadtxt(vis_path+'/lick_index.txt',usecols=[-1],dtype=str,comments='#')

        self.index_list = []
        for i,name in enumerate(names):
            ind = {}
            ind["name"] = name
            ind["type"] = "EW"
            ind["continuum"] = [[lick_index[i,2], lick_index[i,3]], [lick_index[i,4], lick_index[i,5]]]
            ind["feature"] = [lick_index[i,0], lick_index[i,1]]
            ind["unit"] = units[i]
            self.index_list.append(ind)

        self.GUI(figsize=(15,12), index_plot=index_plot)
        
    def ribbon_plot(self, parameter, range=None, log_space=False, nlines=10, lw=1, alpha=1.0, 
                     show=True, figsize=(13,10), cmap='viridis', reverse=False):
        """
        Creates a 2-panel spectrum and residual plot where the varying spectra of a 
        galaxy with the defined initial components under the change of one specified
        parameter. Spectra are in lines that match with the specified colourmap, 
        shown as a colorbar on the side.
        When range is not specified, it defaults to the slider ranges in slider_params.
        Putting log_space=True spaces the parameter values for the spectra in log,
        instead of the default linear spacing.
        Putting reverse=True inverts the zorder and colour scheme of the lines, 
        included for cases where it is easier to spot variations by plotting the 
        spectra of higher parameter values at lower zorder.
        Returns figure and axes (in a list) for further manipulation by the user if needed.
        """
        # create a dummy model
        init_input_logM, init_sfh, init_custom_sfh = utils.create_sfh(self.init_comp)
        model = pipes.model_galaxy(
            utils.make_pipes_components(self.init_comp, init_custom_sfh),
            spec_wavs=self.wavelengths)
        
        # calculate the ticks at which to vary the chosen parameter
        if range is None:
            range = slider_params.slider_lib[parameter]['lims']
        if log_space:
            param_ticks = 10**np.linspace(np.log10(range[0]),np.log10(range[1]),nlines)
        else:
            param_ticks = np.linspace(range[0],range[1],nlines)
        if reverse:
            param_ticks = param_ticks[::-1]
        
        # loop through parameter ticks and update model to get spectrum and residuals to plot
        spectrums = []
        ymaxs = []
        residuals = []
        for i,ticki in enumerate(param_ticks):
            new_init_comp = copy.deepcopy(self.init_comp)
            if ':' in parameter:
                new_init_comp[parameter.split(':')[0]][parameter.split(':')[1]] = ticki
            else:
                new_init_comp[parameter] = ticki
    
            init_input_logM, init_sfh, init_custom_sfh = utils.create_sfh(new_init_comp)
    
            model.update(utils.make_pipes_components(new_init_comp, init_custom_sfh))
            zoom_in_spec = model.spectrum[np.where((model.spectrum[:,0] >= self.init_spec_lim[0]) & 
                                                   (model.spectrum[:,0] <= self.init_spec_lim[1]))]
            spectrums.append(zoom_in_spec.copy())
    
            #figure out the optimal y scale to use, first calculate all ymaxes
            ymaxs.append(1.05*np.max(spectrums[-1][:, 1]))
    
            #calculate residuals (full range)
            run_med = utils.running_median(model.spectrum[:,0], model.spectrum[:,1], width=self.median_width)
            residual = model.spectrum[:,1] / run_med
            zoom_in_res = residual[np.where((model.spectrum[:,0] >= self.init_spec_lim[0]) & 
                                            (model.spectrum[:,0] <= self.init_spec_lim[1]))]
            residuals.append(zoom_in_res.copy())
        
        # extract colours from the chosen colourmap
        if reverse:
            colormap = matplotlib.cm.get_cmap(cmap+'_r')
        else:
            colormap = matplotlib.cm.get_cmap(cmap)
        colors = matplotlib.cm.get_cmap(cmap)((param_ticks-param_ticks[0])/(param_ticks[-1]-param_ticks[0]))
        
        # now create the figure
        fig = plt.figure(figsize=figsize)
    
        gs1 = matplotlib.gridspec.GridSpec(7, 1, hspace=0., wspace=0.)
        ax1 = plt.subplot(gs1[:5])       #main spectrum plot
        ax2 = plt.subplot(gs1[5:])       #residual plot
    
        ymax = max(ymaxs)
        y_scale = int(np.log10(ymax))-1
        res_lims = []
        for i,spectrum in enumerate(spectrums):
            ax1.plot(spectrum[:, 0], spectrum[:, 1]*10**-y_scale,
                    color=colors[i], zorder=4, lw=lw, alpha=alpha)
            ax2.plot(spectrum[:, 0], residuals[i], color=colors[i], lw=lw, zorder=1, alpha=alpha)
            res_span = max(residuals[i]) - min(residuals[i])
            res_lims.append([min(residuals[i])-0.1*res_span, max(residuals[i])+0.1*res_span])
    
        # Sort out spectrum limits and axis labels
        ax1.set_ylim(0., ymax*10**-y_scale)
        ax1.set_xlim(self.init_spec_lim)
        ax1.set_xticks([])
        pipes.plotting.auto_axis_label(ax1, y_scale, z_non_zero=True)
    
        # add residual value guidelines and labels
        ax2.axhline(1, color="black", ls="--", lw=1, zorder=0)
        ax2.axhline(1.5, color="black", ls=":", lw=1, zorder=0)
        ax2.annotate('1.5x', [0.98*(self.init_spec_lim[1]-self.init_spec_lim[0])+self.init_spec_lim[0], 1.5], 
                     ha='center', va='center')
        ax2.axhline(0.5, color="black", ls=":", lw=1, zorder=0)
        ax2.annotate('0.5x', [0.98*(self.init_spec_lim[1]-self.init_spec_lim[0])+self.init_spec_lim[0], 0.5], 
                     ha='center', va='center')
    
        # Sort out residual limits and axis labels
        ax2.set_xlim(self.init_spec_lim)
        pipes.plotting.auto_x_ticks(ax2)
        pipes.plotting.auto_axis_label(ax2, -1, z_non_zero=True)
        ax2.set_ylabel('flux/\ncontinuum')
        # rescale the y axis to be determined only by the residuals in frame
        res_lims = np.array(res_lims)
        ax2.set_ylim([min(res_lims[:,0]), max(res_lims[:,1])])
    
        plt.subplots_adjust(right=0.90)
        cax = plt.axes([0.91, 0.14, 0.02, 0.72])
        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=colormap,
                                                  norm=matplotlib.colors.Normalize(vmin=range[0], 
                                                                                   vmax=range[1])), 
                            cax=cax)
        cbar.set_label(slider_params.slider_lib[parameter]['label'], rotation=270)
    
        if show:
            plt.show()
            
        return fig, [ax1, ax2]
