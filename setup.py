from setuptools import setup
from glob import glob
import re
import io
import os


# Get version strip
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

setup(
    name="BactVarMiner",
    version="1.0.0",
    author="Adam Dinan",
    author_email="adam1989ie@gmail.com",
    description="A command line tool for annotating bacterial genetic variants",
    packages=['BactVarMiner'],
    url="https://github.com/adamd3/BactVarMiner",
    python_requires='>=3.7.0',
    # install_requires=[]
    keywords='bacterial variant genetic annotation burden scoring',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts':
        ['BactVarMiner = BactVarMiner.__main__:main',],
    },
    extras_require={"test": "pytest"},
)
