import os
import sys
import psutil
import platform


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

    else:
        setup_local()


def setup_local():
    """
    Setup the local environment.

    This function adds the path to the momics package to the sys.path.

    Note: I do not install the package via pip install -e, I rather add the path to the package to the sys.path
    -> faster prototyping of the momics package
    """
    if platform.system() == "Linux":
        print("local Linux")
        sys.path.append("/home/davidp/python_projects/marine_omics/marine-omics")
    elif platform.system() == "Windows":
        print("local Windows")
        sys.path.append(
            "C:/Users/David Palecek/Documents/Python_projects/marine_omics/marine-omics"
        )
    else:
        raise NotImplementedError


def setup_ipython():
    """
    Setup the IPython environment.

    This function installs the momics package and other dependencies for the IPython environment.
    """
    if "google.colab" in str(get_ipython()):
        print("Google Colab")

        # !git clone https://github.com/palec87/momics-demos.git
        # this gives access to utils module
        # NOOO, this you need to do in the NB
        # try:
        #     os.system('git clone https://github.com/palec87/momics-demos.git')
        #     print(f"Repository cloned")
        # except OSError as e:
        #     print(f"An error occurred while cloning the repository: {e}")

        # clone and install momics
        try:
            os.system("git clone https://github.com/palec87/marine-omics.git")
            print(f"Repository marine-omics cloned")
        except OSError as e:
            print(f"An error occurred while cloning the repository: {e}")

        try:
            os.system("pip install git+https://github.com/palec87/marine-omics.git")
            print(f"momics installed")
        except OSError as e:
            print(f"An error occurred while installing momics: {e}")

        # !pip install scikit-bio
        # try:
        #     os.system("pip install scikit-bio")
        #     print(f"scikit-bio installed")
        # except OSError as e:
        #     print(f"An error occurred while installing scikit-bio: {e}")

        # !pip install panel hvplot
        try:
            os.system("pip install panel hvplot")
            print(f"panel and hvplot installed")
        except OSError as e:
            print(f"An error occurred while installing panel and hvplot: {e}")

    else:
        # assume local jupyterlab which has all the dependencies installed
        setup_local()


def is_ipython():
    """
    Check if the current environment is IPython.

    Returns:
        bool: True if the current environment is IPython, False otherwise.
    """
    # This is for the case when the script is run from the Jupyter notebook
    if "ipykernel" in sys.modules:
        from IPython import get_ipython

        return True
    else:
        print("Not IPython setup")
        return False


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
