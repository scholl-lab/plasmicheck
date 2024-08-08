from setuptools import setup, find_packages

setup(
    name="plasmicheck",
    version="0.3.0",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pysam",
        "jinja2",
        "weasyprint",
        "matplotlib",
        "seaborn",
        "pandas"
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
