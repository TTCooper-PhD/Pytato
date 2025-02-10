from setuptools import setup, find_packages

setup(
    name="PytatoDIA",
    version="0.5.3",
    author="Tyler T. Cooper,Ph.D.",
    author_email="tcoope2@gmail.com",
    description="A Toolbox for Proteomic Analyses and Data Visualization",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/TTCooper-PhD/Pytato",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "matplotlib>=3.3.4",
        "seaborn>=0.11.1",
        "scipy>=1.5.4",
        "biopython>=1.78",
        "requests>=2.25.1",
        "pandas>=1.2.0",
        "numpy>=1.19.0",
        "logomaker>=0.8",
        "pyteomics>=4.5"
    ],
)
