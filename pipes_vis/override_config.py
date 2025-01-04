import numpy as np
import os
from astropy.io import fits

def override_config(pipes):
    """ 
    Manually changing some of the config from bagpipes to desirable ones that allow
    for the GUI to work properly
    """
    pipes.config.max_redshift = 10
    pipes.config.R_spec = 1800
    
    # Redshift points for the IGM grid.
    pipes.config.igm_redshifts = np.arange(0.0, pipes.config.max_redshift + 0.01, 0.01)
    
    # If the IGM grid has not yet been calculated, calculate it now.
    if not os.path.exists(pipes.config.grid_dir + "/d_igm_grid_inoue14.fits"):
        pipes.config.igm_inoue2014.make_table(pipes.config.igm_redshifts, pipes.config.igm_wavelengths)
    
    else:
        # Check that the wavelengths and redshifts in the igm file are right
        pipes.config.igm_file = fits.open(pipes.config.grid_dir + "/d_igm_grid_inoue14.fits")
    
        if len(pipes.config.igm_file) != 4:
            pipes.config.igm_inoue2014.make_table(pipes.config.igm_redshifts, pipes.config.igm_wavelengths)
    
        else:
            pipes.config.wav_check = np.min(pipes.config.igm_file[2].data == pipes.config.igm_wavelengths)
            pipes.config.z_check = pipes.config.igm_file[3].data.shape[0] == pipes.config.igm_redshifts.shape[0]
    
            if not pipes.config.wav_check or not pipes.config.z_check:
                pipes.config.igm_inoue2014.make_table(pipes.config.igm_redshifts, pipes.config.igm_wavelengths)
    
    # 2D numpy array containing the IGM attenuation grid.
    pipes.config.raw_igm_grid = fits.open(pipes.config.grid_dir + "/d_igm_grid_inoue14.fits")[1].data
