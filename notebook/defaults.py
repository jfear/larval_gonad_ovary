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
from larval_gonad_ovary.notebook import Nb
from larval_gonad_ovary.plotting import make_figs
from larval_gonad_ovary.config import memory

# Setup notebook
nbconfig = Nb.setup_notebook()
