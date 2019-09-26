from setuptools import setup, find_packages
from codecs import open
import os
import re

with open("README.rst", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

def get_version():
    here = os.path.abspath(os.path.dirname(__file__))
    version_file = os.path.join(here, "sage_analysis", "__version__.py")

    with open(version_file, "r") as vf:
        lines = vf.read()
        version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", lines, re.M).group(1)
        return version


sage_analysis_version = get_version()

setup(
    name="sage_analysis",
    version=sage_analysis_version,
    author="Jacob Seiler",
    author_email="jacob.seiler94@gmail.com",
    url="https://github.com/sage-home/sage-analysis",
    download_url="https://pypi.org/project/sage-analysis",
    project_urls={
        "Documentation": "https://sage-model.readthedocs.io",
        "Source Code": "https://github.com/sage-home/sage-analysis"
        },
    description="Ingests, analyses, and plots SAGE data products.",
    long_description=long_description,
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        ],
    package_dir={"sage_analysis": "sage_analysis"},
    packages=find_packages(),
    keywords=("SAGE redshift astronomy astrophysics galaxy semi-analytic model")
)
