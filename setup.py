import os
from setuptools import setup, find_packages

# Load version from version.py
version = {}
with open(os.path.join('plasmicheck', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name="plasmicheck",
    version=version['__version__'],  # Use the version from version.py
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "biopython>=1.84",
        "pysam>=0.22.1",
        "jinja2>=3.0.0",
        "weasyprint>=62.3",
        "matplotlib>=3.9.1.post1",
        "seaborn>=0.13.2",
        "pandas>=2.2.2",
        "scipy>=1.13.1",
        "plotly>=5.23.0",
        "statsmodels>=0.14.2",
        "numpy>=2.0.1",
    ],
    entry_points={
        "console_scripts": [
            "plasmicheck=plasmicheck.cli:main",
        ],
    },
    author="Bernt Popp",
    author_email="bernt.popp.md@gmail.com",
    description="plasmicheck: Detect and quantify plasmid DNA contamination in sequencing data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/berntpopp/plasmicheck",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
        ],
    },
)
