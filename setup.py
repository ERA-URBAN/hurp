import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def conf_path(name):
    if sys.prefix == '/usr':
        conf_path = os.path.join('/etc', name)
    else:
        conf_path = os.path.join(sys.prefix, 'etc', name)
    return conf_path

setup(
    name = "hurp",
    version = "0.0.1",
    author = "Ronald van Haren, Michiel van der Molen",
    author_email = "r.vanharen@esciencecenter.nl",
    description = ("A python library to create WRFCHEM input files for WRF"),
    license = "Apache 2.0",
    keywords = "WRF CHEM",
    url = "https://github.com/ERA-URBAN/hurp",
    packages=['hurp'],
    scripts=['hurp/scripts/hurp'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved ::Apache Software License",
    ],
)
