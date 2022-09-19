import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pipes-vis",
    version="0.4.2",
    author="Ho-Hin Leung",
    author_email="hhl1@st-andrews.ac.uk",
    description="Small interactive GUI/visualizer tool for SPS spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/HinLeung622/pipes_vis",
    project_urls={
        "Bug Tracker": "https://github.com/HinLeung622/pipes_vis/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    
    install_requires=["bagpipes>=0.8.6", "numpy>=1.14.2",
                      "astropy", "matplotlib>=2.2.2", "scipy"],
)
