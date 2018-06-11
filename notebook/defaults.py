# Imports
import os
import sys
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from lib.notebook import Nb
from lib.plotting import make_figs
from lib.config import memory

# Setup notebook
nbconfig = Nb.setup_notebook()
