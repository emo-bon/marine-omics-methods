import os
import sys
import psutil
import logging
from IPython import get_ipython


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
