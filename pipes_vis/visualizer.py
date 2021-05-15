import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox

try:
    import bagpipes as pipes
except ImportError:
    print('BAGPIPES not installed')
    #from . import pipes_source as pipes

from . import utils
from . import plotting
from . import slider_params


class visualizer:
    """
    Object that holds and does the under-the-hood work for showing the galaxy SFH/parameters
    -> spectrum GUI visualizer.
    """
    def __init__(self, init_components):
        self.init_comp = init_components
        self.spec_lim = self.init_comp['spec_lim'].copy()
        # full spectrum wavelengths
        wavelengths_fine = np.linspace(1000,10000, 10000)
        wavelengths_coarse = np.linspace(10000,50000,2000)[1:]
        self.wavelengths = np.concatenate([wavelengths_fine, wavelengths_coarse])
        # colours used for the various plot lines
        self.plot_colors = {'sfh':"black", 'z':"red", 'spectrum':"sandybrown", 
                            'continuum':"black", 'zoom':"black"}
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
        self.sub_ax_arg = [.69, .25, .2, .11]          # location and scaling of full spectrum subplot
        self.bottom_adjust_val = 0.36                  # scaling value used on all plots to make space
        self.top_adjust_val = 0.95                     # scaling value used on all plots to make space
                                                       # for the sliders
        self.slider_colors = {'sfh':['lightsteelblue','bisque','honeydew','mistyrose','thistle'],
                              'misc':'lightgoldenrodyellow'}      # colour of the sliders
        self.reset_button_arg = [0.8, 0.025, 0.1, 0.04]# location and scaling of the reset button
        
    def static_plot(self, show=True, figsize=(13,9)):
        """ 
        Creates the figure, lines, texts and annotations. Returns figure and axes (in a list)
        for further manipulation by the user if needed.
        """
        self.fig = plt.figure(figsize=figsize)

        gs1 = matplotlib.gridspec.GridSpec(13, 1, hspace=0., wspace=0.)
        self.ax1 = plt.subplot(gs1[0:5])        #SFH plot
        self.ax2 = plt.subplot(gs1[6:11])       #main spectrum plot
        self.ax3 = plt.subplot(gs1[11:])        #residual plot

        init_input_logM, init_sfh, init_custom_sfh = utils.create_sfh(self.init_comp)

        self.sfh_line, self.z_line, self.z_text, self.input_logM_text, self.bad_sfh_text \
            = plotting.add_sfh_plot(self.init_comp, init_sfh, init_input_logM, self.ax1,
                                    sfh_color=self.plot_colors['sfh'], z_line_color=self.plot_colors['z'])

        self.model = pipes.model_galaxy(utils.make_pipes_components(self.init_comp, init_input_logM, init_custom_sfh),
                                        spec_wavs=self.wavelengths)

        # full spectrum in inset
        self.sub_ax = plt.axes(self.sub_ax_arg)
        sub_y_scale_spec,self.sub_spec_line = plotting.add_bp_spectrum(self.model.spectrum, self.sub_ax, sub=True,
                                                                       color=self.plot_colors['spectrum'])
        self.spec_zoom_poly = self.sub_ax.fill_between(self.spec_lim, [0]*2, [20]*2, color=self.plot_colors['zoom'], 
                                                       alpha=0.1)

        # the main spectrum plot
        self.spec_line, self.run_med_line, self.overflow_text, y_scale_spec \
            = plotting.add_main_spec(self.model, self.ax2, self.spec_lim, median_width=self.median_width, 
                                     color=self.plot_colors['spectrum'],
                                     continuum_color=self.plot_colors['continuum'])

        # the residual plot
        self.res_line = plotting.add_residual(self.model, self.ax3, self.spec_lim, median_width=self.median_width, 
                                              color=self.spec_line.get_color())
        
        if show:
            plt.show()
            
        return self.fig, [self.ax1, self.ax2, self.ax3, self.sub_ax]
        
    def GUI(self, figsize=(13,12)):
        """ 
        Creates the figure, lines, texts and annotations that will be manipulated, and also the 
        interactive elements.
        """
        self.static_plot(show=False, figsize=figsize)

        # adjust the main plots to make room for the sliders
        plt.subplots_adjust(bottom=self.bottom_adjust_val, top=self.top_adjust_val)
        default_top_margin = 1-plt.rcParams["figure.subplot.top"]
        default_bottom_margin = plt.rcParams["figure.subplot.bottom"]
        new_top_margin = 1-self.top_adjust_val
        new_bottom_margin = self.bottom_adjust_val
        ratio = (1-new_top_margin-new_bottom_margin) / (1-default_top_margin-default_bottom_margin)
        
        l, b, w, h = self.sub_ax.get_position().bounds
        new_b = 1-new_top_margin-(1-b-default_top_margin)*ratio
        new_h = h*ratio
        self.sub_ax.set_position([l, new_b, w, new_h])

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
        self.spec_min_box = TextBox(self.ax_spec_min, r'$\lambda_{min}$', initial=self.spec_lim[0])
        self.spec_min_box.on_submit(self.submit_min)

        self.ax_spec_max = plt.axes([current_right_column_pos[0]+textbox_width+self.textbox_gap, 
                                     current_right_column_pos[1], textbox_width, self.height])
        self.spec_max_box = TextBox(self.ax_spec_max, r'$\lambda_{max}$', initial=self.spec_lim[1])
        self.spec_max_box.on_submit(self.submit_max)

        # register the update function with each slider
        for key in self.sliders.keys():
            self.sliders[key].on_changed(self.update)

        # Create a `matplotlib.widgets.Button` to reset all sliders to initial values.
        self.resetax = plt.axes(self.reset_button_arg)
        self.reset_button = Button(self.resetax, 'Reset', color=self.slider_colors['misc'], hovercolor='0.975')
        
        self.reset_button.on_clicked(self.reset)

        plt.show()
        
    def make_one_slider(self, x_pos, y_pos, width, height, label, lims, init_val, bg_color):
        """ makes a single slider """
        ax_slider = plt.axes([x_pos, y_pos, width, height], facecolor=bg_color)
        the_slider = Slider(
            ax=ax_slider,
            label=label,
            valmin=lims[0],
            valmax=lims[1],
            valinit=init_val,
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
        custom_extremes = []
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
                        if key in sfh_types:
                            required_components.append(key+':'+sfh_key)
                        elif key[:-1] in sfh_types:
                            required_components.append(key[:-1]+':'+sfh_key)
                        slider_names.append(key+':'+sfh_key)
                        init_vals.append(sfh_dict[sfh_key])
                        try:
                            custom_extremes.append(sfh_dict[sfh_key+'_lims'])
                        except KeyError:
                            custom_extremes.append(None)
                            
            elif key in ['dust', 'nebular']:
                for sub_key in self.init_comp[key].keys():
                    if sub_key[-5:] != '_lims':
                        if sub_key != 'type':
                            required_components.append(key+':'+sub_key)
                            slider_names.append(key+':'+sub_key)
                            init_vals.append(self.init_comp[key][sub_key])
                            try:
                                custom_extremes.append(self.init_comp[key][sub_key+'_lims'])
                            except KeyError:
                                custom_extremes.append(None)
                                
            elif key != 'spec_lim' and key[-5:] != '_lims':
                required_components.append(key)
                slider_names.append(key)
                init_vals.append(self.init_comp[key])
                try:
                    custom_extremes.append(self.init_comp[key+'_lims'])
                except KeyError:
                    custom_extremes.append(None)
            
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
            
            if custom_extremes[i] is not None:
                c_dict['lims'] = custom_extremes[i]
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
    
    def update_spec(self, spectrum, ax, spec_line, sub=False, overflow_text=None, change_xlims=False):
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

    def update(self, val):
        """ The function to be called anytime a slider's value changes """
        #create updated components dictionary
        sfh_dict_list = []
        sfh_types = dir(pipes.models.star_formation_history)
        for key in self.init_comp.keys():
            if key in sfh_types or key[:-1] in sfh_types:
                sfh_dict = {}
                for sfh_key in self.init_comp[key].keys():
                    if sfh_key[-5:] != '_lims':
                        sfh_dict[sfh_key] = self.sliders[key+':'+sfh_key].val
                sfh_dict_list.append(sfh_dict)
        
        other_modules = {}
        for module in ['dust', 'nebular']:
            if module in self.init_comp.keys():
                module_dict = {}
                for key in self.init_comp[module].keys():
                    if key == 'type':
                        module_dict["type"] = self.init_comp[module]['type']
                    elif key[-5:] != '_lims':
                        module_dict[key] = self.sliders[module+':'+key].val
                other_modules[module] = module_dict

        self.new_comp = {}
        sfh_index = 0
        for key in self.init_comp.keys():
            if key in sfh_types or key[:-1] in sfh_types:
                self.new_comp[key] = sfh_dict_list[sfh_index]
                sfh_index += 1
            elif key == 'dust':
                self.new_comp["dust"] = other_modules['dust']
            elif key == 'nebular':
                self.new_comp["nebular"] = other_modules['nebular']
            elif key != 'spec_lim' and key[-5:] != '_lims':
                self.new_comp[key] = self.sliders[key].val

        #update sfh plot
        age_at_z = utils.cosmo.age(self.new_comp["redshift"]).value
        self.z_line.set_xdata([age_at_z,age_at_z])
        input_logM, vary_sfh, custom_sfh = utils.create_sfh(self.new_comp)
        self.sfh_line.set_ydata(vary_sfh[1])
        if max(vary_sfh[1])>0:
            self.ax1.set_ylim(top=1.05*max(vary_sfh[1]))
            self.bad_sfh_text.set_alpha(0.0)
        else:
            self.ax1.set_ylim(top=1)
            self.bad_sfh_text.set_alpha(1.0)
        self.z_text.remove()
        self.z_text = self.ax1.annotate(str(np.round(self.new_comp["redshift"],3)), 
                                        [age_at_z, 0.92*self.ax1.get_ylim()[1]], color=self.z_line.get_color())
        self.input_logM_text.set_text('pre-obs log10M = '+str(np.round(input_logM,2)))

        #update model components
        self.model.update(utils.make_pipes_components(self.new_comp, input_logM, custom_sfh))
        
        run_med = utils.running_median(self.model.spectrum[:,0], self.model.spectrum[:,1], 
                                       width=self.median_width)
        residual = self.model.spectrum[:,1] / run_med

        #update subplot full spectrum
        self.update_spec(self.model.spectrum, self.sub_ax, self.sub_spec_line, sub=True)

        #update main spectrum plot
        zoom_in_spec = self.model.spectrum[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                                    (self.model.spectrum[:,0] <= self.spec_lim[1]))]
        y_scale_spec = self.update_spec(zoom_in_spec, self.ax2, self.spec_line, overflow_text=self.overflow_text)
        if y_scale_spec is not None:
            self.run_med_line.set_ydata(run_med*10**-y_scale_spec)

            #update residual plot
            self.res_line.set_ydata(residual)
            in_range_res = residual[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                             (self.model.spectrum[:,0] <= self.spec_lim[1]))]
            res_span = max(in_range_res) - min(in_range_res)
            self.ax3.set_ylim([min(in_range_res)-0.1*res_span, max(in_range_res)+0.1*res_span])

        self.fig.canvas.draw_idle()
        

    def submit_min(self, text):
        """
        adjust the main spectrum plot's limits and the full spec's shaded region
        on submit of a new wavelength lower limit value
        """
        self.spec_lim[0] = eval(text)
        zoom_in_spec = self.model.spectrum[np.where((self.model.spectrum[:,0] >= self.spec_lim[0]) & 
                                               (self.model.spectrum[:,0] <= self.spec_lim[1]))]
        self.update_spec(zoom_in_spec, self.ax2, self.spec_line, change_xlims=True)
        self.ax3.set_xlim(self.spec_lim)
        pipes.plotting.auto_x_ticks(self.ax3)
        self.spec_zoom_poly.set_verts([[[self.spec_lim[0],0],
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
        self.update_spec(zoom_in_spec, self.ax2, self.spec_line, change_xlims=True)
        self.ax3.set_xlim(self.spec_lim)
        pipes.plotting.auto_x_ticks(self.ax3)
        self.spec_zoom_poly.set_verts([[[self.spec_lim[0],0],
                                   [self.spec_lim[0],20],
                                   [self.spec_lim[1],20],
                                   [self.spec_lim[1],0]]])
        self.fig.canvas.draw_idle()
        
    def reset(self, event):
        """ resets everything to the default states """
        for key in self.sliders.keys():
            self.sliders[key].reset()
        self.spec_min_box.set_val(self.init_comp['spec_lim'][0])
        self.spec_max_box.set_val(self.init_comp['spec_lim'][1])
