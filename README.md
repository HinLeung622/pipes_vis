# pipes_vis
A small, interactive GUI/visualizer tool for SPS spectra, powered by Adam Carnall's Bayesian Analysis of Galaxies for Physical Inference and Parameter EStimation [(BAGPIPES)](https://bagpipes.readthedocs.io/en/latest/) <br>

A ["Mini-paper"](https://iopscience.iop.org/article/10.3847/2515-5172/ac1680) about pipes_vis is published in the Research Notes of the AAS. Please cite this paper if you use pipes_vis in your publications.<br>

To install, type `pip install pipes_vis`. This package depends on Bagpipes, which has non-python packages like multinest as dependencies. If you only wish to use pipes_vis, you need not go through the process to get all non-python dependencies of bagpipes as they are only used for SED fitting rather than modelling. Just ignore complains from Bagpipes when importing. However, if you wish to use Bagpipes for SED fitting, follow the tutorial in Bagpipes' docs to get those non-python packages.<br>

You can find examples and documentations of the tool in [this jupyter notebook](https://github.com/HinLeung622/pipes_vis/blob/main/pipes_vis_example.ipynb)

Dependencies: bagpipes, numpy, matplotlib, astropy, scipy (visualizer still works even if you do not install multinest and pymultinest as instructed by bagpipes' installation) <br>

For requests regarding making the visualizer avaliable as a web interactive tool much like [sengi](https://www.christopherlovell.co.uk/sengi/), it will likely require rewriting large portions of the code in javascript or some other web-purposed languages. Please contact me if you are interested in doing this/want this to be done.<br>

![alt text](https://github.com/HinLeung622/pipes_vis/blob/main/readme_image.png?raw=true)
