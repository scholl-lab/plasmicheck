import json
from setuptools import setup, find_packages

# Load configuration from JSON file
with open('config.json', 'r') as config_file:
    config = json.load(config_file)

setup(
    name="plasmicheck",
    version=config['version'],
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pysam",
        "jinja2",
        "weasyprint",
        "matplotlib",
        "seaborn",
        "pandas",
        "scipy",
        "plotly",
        "statsmodels"
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
)
