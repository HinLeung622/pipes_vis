from . import utils
from . import plotting
from . import slider_params
from . import override_config
from . import fix_wavelength_range_bug
from .visualizer import visualizer
print(
      'If bagpipes is complaining about a missing installation of multinest and/or pymultinest, '+
      'it means you will not be able to use the fitting half of bagpipes. But that will not '+
      'affect the usage of pipes_vis. However, if you do intend to use bagpipes for SED fitting, '+
      'do follow the relevent instructions to install the required packages.'
      )
