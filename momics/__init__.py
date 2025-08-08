import importlib as _importlib
import logging

# from momics.version import version as __version__

# This goes into your library somewhere
logging.getLogger("momics").addHandler(logging.NullHandler())

submoduless = [
    "diversity",
    "galaxy",
    "loader",
    "metadata",
    "networks",
    "plotting",
    "stats",
    "taxonomy",
]
