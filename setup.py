#! /usr/bin/env python
#
# 2022 Vincenzo Marra

DESCRIPTION = "NeurophysioTools: tools for Neurophysiological signals in Python"
DISTNAME = "neurophysiotools"
MAINTAINER = "Vincenzo Marra"
MAINTAINER_EMAIL = "raphaelvallat9@gmail.com"
URL = "https://github.com/enzomarra/neurophysiotools"
LICENSE = "MIT license"
DOWNLOAD_URL = "https://github.com/enzomarra/neurophysiotools"
VERSION = "0.0.1"
PACKAGE_DATA = {}

try:
    from setuptools import setup

    _has_setuptools = True
except ImportError:
    from distutils.core import setup


def check_dependencies():
    install_requires = []

    try:
        import numpy
    except ImportError:
        install_requires.append("numpy")
    try:
        import scipy
    except ImportError:
        install_requires.append("scipy")


    return install_requires


if __name__ == "__main__":

    install_requires = check_dependencies()

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        install_requires=install_requires,
        include_package_data=True,
        packages=["neurophysiotools"],
        package_data=PACKAGE_DATA,
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python",
            "License :: OSI Approved :: MIT License",
            "Topic :: Scientific/Engineering",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
    )
