import bagpipes as pipes
from bagpipes import config
import numpy as np

# fix the bug that is fixed in bagpipes commit e8fc99d but not yet publically available

def _get_wavelength_sampling(self):
    """ Calculate the optimal wavelength sampling for the model
    given the required resolution values specified in the config
    file. The way this is done is key to the speed of the code. """
    
    max_z = config.max_redshift
    
    if self.spec_wavs is None:
        self.max_wavs = [(self.filter_set.min_phot_wav
                          / (1.+max_z)),
                         1.01*self.filter_set.max_phot_wav, 10**8]
    
        self.R = [config.R_other, config.R_phot, config.R_other]
    
    elif self.filt_list is None:
        self.max_wavs = [self.spec_wavs[0]/(1.+max_z),
                         self.spec_wavs[-1], 10**8]
    
        self.R = [config.R_other, config.R_spec, config.R_other]
    
    else:
        if (self.spec_wavs[0] > self.filter_set.min_phot_wav
                and self.spec_wavs[-1] < self.filter_set.max_phot_wav):
    
            self.max_wavs = [self.filter_set.min_phot_wav/(1.+max_z),
                             self.spec_wavs[0]/(1.+max_z),
                             self.spec_wavs[-1],
                             self.filter_set.max_phot_wav, 10**8]
    
            self.R = [config.R_other, config.R_phot, config.R_spec,
                      config.R_phot, config.R_other]
    
        elif (self.spec_wavs[0] < self.filter_set.min_phot_wav
              and self.spec_wavs[-1] < self.filter_set.max_phot_wav):
    
            self.max_wavs = [self.spec_wavs[0]/(1.+max_z),
                             self.spec_wavs[-1],
                             self.filter_set.max_phot_wav, 10**8]
    
            self.R = [config.R_other, config.R_spec,
                      config.R_phot, config.R_other]
    
        elif (self.spec_wavs[0] > self.filter_set.min_phot_wav
              and self.spec_wavs[-1] > self.filter_set.max_phot_wav):
    
            self.max_wavs = [self.filter_set.min_phot_wav/(1.+max_z),
                             self.spec_wavs[0]/(1.+max_z),
                             self.spec_wavs[-1], 10**8]
    
            self.R = [config.R_other, config.R_phot,
                      config.R_spec, config.R_other]
    
        elif (self.spec_wavs[0] < self.filter_set.min_phot_wav
              and self.spec_wavs[-1] > self.filter_set.max_phot_wav):
    
            self.max_wavs = [self.spec_wavs[0]/(1.+max_z),
                             self.spec_wavs[-1], 10**8]
    
            self.R = [config.R_other, config.R_spec, config.R_other]
    
    # Generate the desired wavelength sampling.
    x = [1.]
    
    for i in range(len(self.R)):
        if i == len(self.R)-1 or self.R[i] > self.R[i+1]:
            while x[-1] < self.max_wavs[i]:
                x.append(x[-1]*(1.+0.5/self.R[i]))
    
        else:
            while x[-1]*(1.+0.5/self.R[i]) < self.max_wavs[i]:
                x.append(x[-1]*(1.+0.5/self.R[i]))
    
    return np.array(x)

pipes.models.model_galaxy._get_wavelength_sampling = _get_wavelength_sampling
