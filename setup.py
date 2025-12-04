#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path

print("If you are using ptt to run calculations through MAGEMin, please run 'ptt.install_MAGEMinCalc()' or 'ptt.update_MAGEMinCalc()' to update the associated julia code.")
print("If you notice any bugs/issues please contact me at gleesonm@berkeley.edu")

this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'src', 'petthermotools', '_version.py'), encoding='utf-8') as f:
    exec(f.read())

with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name="petthermotools",
    version=__version__,
    author="Matthew Gleeson",
    author_email="gleesonm@berkeley.edu",
    description="petthermotools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gleesonm1/PetThermoTools",
    package_dir={'': 'src'},  # Optional
    packages=find_packages(where='src', include=['petthermotools', 'petthermotools.*']),  # Required

    package_data={
        # Include all pickle files
        "": ["*.jl"],
    },
    install_requires=[
            'pandas',
            'numpy',
            'matplotlib',
            'scikit-learn',
            'scipy',
            'tinynumpy',
            'shapely',
            'tqdm'
            ],

    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)
