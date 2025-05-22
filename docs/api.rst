******
API
******

.. module:: momics

This is the API documentation for the momics package

Loading module
==================
The loading module provides utilities for importing and handling various data formats used in marine omics workflows.
Metadata are handled in separate module `momics.metadata`.

Parquet data
-----------------
This submodule contains functions for loading and processing data stored in Parquet format.

.. automodule:: momics.loader.parquets
    :members:
    :show-inheritance:

Ro-crates
-----------------
This submodule provides tools for working with RO-Crate metadata packages.

.. automodule:: momics.loader.ro_crates
    :members:
    :show-inheritance:

Utils
-----------------
This submodule includes helper functions used during the data loading process.

.. automodule:: momics.loader.utils
    :members:
    :show-inheritance:


Diversity module
==================
This module offers methods for calculating and analyzing biodiversity metrics from omics data.

.. automodule:: momics.diversity
    :members:
    :show-inheritance:

Panel dashboard module
========================
This module provides some of the utilities for building interactive dashboards using the Panel library. These are Specific
for FAIR-EASE use case and more methods and widgets can be found directly in the demo workflow notebooks [here](https://github.com/emo-bon/momics-demos).

.. automodule:: momics.panel_utils
    :members:
    :show-inheritance:

Plotting module
==================
This module contains functions for visualizing omics data using various plotting libraries, such as `mpl`, `seaborn`, and `hvplot`.

.. automodule:: momics.plotting
    :members:
    :show-inheritance:

Taxonomy module
==================
This module provides tools for handling and analyzing taxonomic information in omics datasets.

.. automodule:: momics.taxonomy
    :members:
    :show-inheritance:

Galaxy integration
==================
This module enables integration with the Galaxy platform for workflow automation and reproducibility. We keep it minimalistic, because we expect to use direct install of Galaxy once demos are deployed in the final VRE.

.. automodule:: momics.galaxy
    :members:
    :show-inheritance:

Utilities of all sorts
==================
This module contains miscellaneous utility functions used throughout the momics package.

.. automodule:: momics.utils
    :members:
    :show-inheritance: