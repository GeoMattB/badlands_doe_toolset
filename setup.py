##############################

##############################

from setuptools import setup, find_packages

import glob
import subprocess
from os import path
import io

sys_includes = []

if __name__ == "__main__":
    setup(
        name="badlands_doe_toolset",
        author="Matt Boyd",
        author_email="mboyd@sydney.edu.au",
        url="https://github.com/GeoMattB",
        version="0.1",
        description="badlands_doe_toolset is a set of tools to help build and analyse Design of Experiment configurations for badlands modelling",
        ext_modules = [],
        packages=['badlands_doe_toolset'],
        package_data={'badlands_doe_toolset':sys_includes},
        data_files=[('badlands_doe_toolset',sys_includes)],
        include_package_data=True,
        install_requires=[
            "badlands>=2.2.3",
            "badlands-companion",
            "similaritymeasures>=0.4.4",
            "doegen",
            "numpy",
            "pandas",
            "networkx",
            "scipy"
        ],
        python_requires=">=3.8",
        classifiers=[
            "Programming Language :: Python :: 3.8",
        ],
    )
