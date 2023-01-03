from setuptools import setup
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
    name="VarPredict",
    version=find_version("VarPredict/__init__.py"),
    author="Adam Dinan",
    author_email="ad866@cam.ac.uk",
    description="A command line tool for modelling gene expression using genetic variants and covariates",
    url="https://github.com/adamd3/VarPredict",
    packages=['VarPredict'],
    python_requires='>=3.7.0',
    install_requires=['matplotlib', 'numpy>=1.21.0', 'pandas>=1.4.0', 'sklearn'],
    keywords='genetic variant expression machine learning model',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
    entry_points={
        'console_scripts':
        ['VarPredict = VarPredict.__main__:main',],
    },
    extras_require={"test": "pytest"},
)
