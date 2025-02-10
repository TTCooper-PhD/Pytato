"""PytatoDIA: A Python package for proteomics data analysis."""

# Import main components
from .core.silico import Silico
from .core.scales import Scales
from .analysis.peptide import Peptide
from .analysis.protein import Protein
from .viz.plots import Plots

# Version info
__version__ = "0.5.3"  # or whatever your current version is