import setuptools
import re
import io
import os


# Get version 
def read(*names, **kwargs):
    with io.open(os.path.join(os.path.dirname(__file__), *names),
                 encoding=kwargs.get("encoding", "utf8")) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setuptools.setup(
    name="BactVarMiner",
    version=find_version("BactVarMiner/__init__.py"),
    author="Adam Dinan",
    author_email="ad866@cam.ac.uk",
    description="A command line tool for modelling gene expression using genetic variants and covariates",
    url="https://github.com/adamd3/BactVarMiner",
    packages=setuptools.find_packages(),
    package_dir={'BactVarMiner':'BactVarMiner'},
    package_data={'BactVarMiner': ['BactVarMiner/*.txt']},
    include_package_data=True,
    python_requires='>=3.7.0',
    install_requires=['ConfigArgParse', 'matplotlib', 'numpy>=1.21.0', 'pandas>=1.4.0', 'sklearn'],
    scripts=['BactVarMiner/BactVarMiner'],
    keywords='variant expression machine learning model',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
