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
    name="DeepBactVAT",
    version="1.0.0",
    author="Adam Dinan, Aaron Weimann",
    author_email="adam1989ie@gmail.com",
    description="A command line tool that uses deep learning to score bacterial genetic variants",
    packages=['DeepBactVAT'],
    url="https://github.com/adamd3/DeepBactVAT",
    python_requires='>=3.8.0',
    # install_requires=[]
    keywords='bacterial variant association eQTL',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts':
        ['DeepBactVAT = DeepBactVAT.__main__:main',],
    },
    extras_require={"test": "pytest"},
)
