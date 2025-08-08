import os
import sys
import psutil
import logging
import pandas as pd
from typing import Tuple
from IPython import get_ipython

import pandas as pd
from momics.taxonomy import (
    prevalence_cutoff,
    remove_high_taxa,
    prevalence_cutoff_taxonomy,
    fill_taxonomy_placeholders,
    pivot_taxonomic_data,
)
from momics.loader import load_parquets_udal
from momics.metadata import (
    get_metadata_udal,
    enhance_metadata,
    clean_metadata,
    merge_source_mat_id_to_data,
)
from momics.constants import COL_NAMES_HASH_EMO_BON_VRE as COL_NAMES_HASH

# logger setup
FORMAT = "%(levelname)s | %(name)s | %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


#####################
# Environment setup #
#####################
# TODO: there needs to be a yaml file to set up a folder structure, hardcoding here is not good :)
# Question: Or should this be part of the momics package?
def init_setup():
    """
    Initializes the setup environment.

    This function checks if the current environment is IPython (such as Google Colab).
    If it is, it runs the setup for IPython environments. Otherwise, it runs the setup
    for local environments.
    """
    # First solve if IPython
    if is_ipython():
        ## For running at GColab, the easiest is to clone and then pip install some deps
        setup_ipython()


def install_colab_packages():
    try:
        os.system("pip install panel hvplot")
        print(f"panel and hvplot installed")
    except OSError as e:
        print(f"An error occurred while installing panel and hvplot: {e}")


def setup_ipython():
    """
    Setup the IPython environment.

    This function installs the momics package and other dependencies for the IPython environment.
    """
    if "google.colab" in str(get_ipython()):
        print("Google Colab")

        # Install ngrok for hosting the dashboard
        try:
            os.system("pip install pyngrok --quiet")
            print("ngrok installed")
        except OSError as e:
            print(f"An error occurred while installing ngrok: {e}")

        # Install the momics package
        install_colab_packages()


def is_ipython():
    # This is for the case when the script is run from the Jupyter notebook
    if "ipykernel" not in sys.modules:
        print("Not an IPython setup")
        return False

    return True


def get_notebook_environment():
    """
    Determine if the notebook is running in VS Code or JupyterLab.

    Returns:
        str: The environment in which the notebook is running ('vscode', 'jupyter:binder', 'jupyter:local' or 'unknown').
    """
    # Check for VS Code environment variable
    if "VSCODE_PID" in os.environ:
        return "vscode"

    elif "JPY_SESSION_NAME" in os.environ:
        if psutil.users() == []:
            print("Binder")
            return "jupyter:binder"
        else:
            print("Local Jupyter")
            return "jupyter:local"
    else:
        return "unknown"


###########
# logging #
###########
FORMAT = "%(levelname)s | %(name)s | %(message)s"  # for logger


def reconfig_logger(format=FORMAT, level=logging.INFO):
    """(Re-)configure logging"""
    logging.basicConfig(format=format, level=level, force=True)

    # removing tarnado access logs
    hn = logging.NullHandler()
    logging.getLogger("tornado.access").addHandler(hn)
    logging.getLogger("tornado.access").propagate = False

    logging.info("Logging.basicConfig completed successfully")


##############
# Dataframes #
##############
# check if two dataframes have the same index name
def check_index_names(df1: pd.DataFrame, df2: pd.DataFrame) -> bool:
    """
    Check if two DataFrames have the same index name.

    Args:
        df1 (pd.DataFrame): The first DataFrame.
        df2 (pd.DataFrame): The second DataFrame.

    Returns:
        bool: True if both DataFrames have the same index name, False otherwise.
    """
    return df1.index.name == df2.index.name


#####################
# Memory management #
#####################
def memory_load():
    """
    Get the memory usage of the current process.

    Returns:
        tuple: A tuple containing:
            - used_gb (float): The amount of memory currently used by the process in gigabytes.
            - total_gb (float): The total amount of memory available in gigabytes.
    """
    used_gb, total_gb = (
        psutil.virtual_memory()[3] / 1000000000,
        psutil.virtual_memory()[0] / 1000000000,
    )
    return used_gb, total_gb


def memory_usage():
    """
    Get the memory usage of the current process.

    Returns:
        list: A list of tuples containing the names of the objects in the current environment
            and their corresponding sizes in bytes.
    """
    # These are the usual ipython objects, including this one you are creating
    ipython_vars = ["In", "Out", "exit", "quit", "get_ipython", "ipython_vars"]

    # Get a sorted list of the objects and their sizes
    mem_list = sorted(
        [
            (x, sys.getsizeof(globals().get(x)))
            for x in dir()
            if not x.startswith("_") and x not in sys.modules and x not in ipython_vars
        ],
        key=lambda x: x[1],
        reverse=True,
    )

    return mem_list


########################
# High-level functions #
########################


def load_and_clean(
    valid_samples: pd.DataFrame = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Load metadata
    full_metadata = get_metadata_udal()

    # filter the metadata only for valid 181 samples
    full_metadata, added_columns = enhance_metadata(full_metadata, valid_samples)

    # LOADing data
    mgf_parquet_dfs = load_parquets_udal()
    mgf_parquet_dfs = merge_source_mat_id_to_data(mgf_parquet_dfs, full_metadata)

    # convert added_columns to a dictionary
    added_columns = {col: col.replace("_", " ") for col in added_columns}

    # extend COL_NAMES_HASH with added columns
    COL_NAMES_HASH.update(added_columns)

    # clean the metadata, rename and remove certain columns for the VRE purposes
    full_metadata = clean_metadata(full_metadata, COL_NAMES_HASH)

    return full_metadata, mgf_parquet_dfs


def taxonomy_common_preprocess01(
    df, high_taxon, mapping, prevalence_cutoff_value, taxonomy_ranks, pivot=False
):
    df1 = fill_taxonomy_placeholders(df, taxonomy_ranks)

    logger.info("Preprocessing taxonomy...")
    if high_taxon != "None":
        bef = df1.shape[0]
        df1 = remove_high_taxa(
            df1, taxonomy_ranks, tax_level=high_taxon, strict=mapping
        )
        aft = df1.shape[0]
        logger.info(f"Removed {bef - aft} high taxa at level: {high_taxon}")

    # low prevalence cutoff
    if pivot:
        df1 = pivot_taxonomic_data(df1)

        # low prevalence cutoff
        df1 = prevalence_cutoff(df1, percent=prevalence_cutoff_value, skip_columns=0)
        return df1
    else:
        return prevalence_cutoff_taxonomy(df1, percent=prevalence_cutoff_value)
