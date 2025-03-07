import os
import logging
from datetime import datetime
from typing import List, Tuple

from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.config import ConfigClient


class RemGalaxy:
    def __init__(self, url_var_name: str, api_key_var_name: str):
        """Initializes the BlueCloud instance.

        Args:
            url_var_name (str): The name of the environment variable containing the Galaxy URL.
            api_key_var_name (str): The name of the environment variable containing the Galaxy API key.

        Raises:
            TypeError: If the provided Galaxy URL or API key is invalid.
            Exception: For any other errors that occur during initialization.
        """
        self.logger = logging.getLogger("BCGalaxy")
        self.logger.setLevel(logging.DEBUG)
        # logging.debug("test")

        # url is a name of the variable saved in the .env file
        # the same for the api_key
        try:
            self.set_galaxy_env(url_var_name, api_key_var_name)
            self.cfg = ConfigClient(self.gi)
            self.logger.info(f"User: {self.cfg.whoami()}")
            self.logger.info(f"Galaxy version: {self.cfg.get_version()}")
        except TypeError:
            self.logger.error("Please provide valid Galaxy URL and API key!")
        except Exception as e:
            self.logger.error(f"Error: {e}")

    def set_galaxy_env(self, url_var_name: str, api_key_var_name: str) -> List:
        """Sets the Galaxy environment variables.

        Args:
            url_var_name (str): The name of the environment variable containing the Galaxy URL.
            api_key_var_name (str): The name of the environment variable containing the Galaxy API key.

        Returns:
            List: A list of environment variables.
        """
        self.url = os.getenv(url_var_name)
        self.api_key = os.getenv(api_key_var_name)
        self.gi = GalaxyInstance(self.url, key=self.api_key)

    def get_datasets_by_key(
        self, key: str, value: str | List[str]
    ) -> Tuple[List, List, List]:
        """Retrieves datasets by a specific key and value.

        Args:
            key (str): The key to filter datasets.
            value (str | List[str]): The value to filter datasets.

        Returns:
            Tuple[List[str], List[str], List[str]]:
                A tuple containing:
                - A list of datasets (tuples) that match the key and value.
                - A list of dataset names that match the key and value.
                - A list of dataset IDs that match the key and value.
        """

        if isinstance(value, list):
            lst_dict = [
                k
                for k in self.gi.datasets.get_datasets()
                if key in k and k[key] in value
            ]
        elif isinstance(value, str):
            lst_dict = [
                k
                for k in self.gi.datasets.get_datasets()
                if key in k and k[key] == value
            ]

        else:
            raise ValueError("Value must be a string or a list of strings.")

        self.ds = [(k["name"], k["id"]) for k in lst_dict]
        self.ds_names = [k[0] for k in self.ds]
        self.ds_ids = [k[1] for k in self.ds]

        return self.ds, self.ds_names, self.ds_ids

    def get_datasets(self, name: str = None):
        """Retrieves datasets with an optional name filter.

        Args:
            name (str, optional): The name to filter datasets. Defaults to None.

        Returns:
            List[str]: A list of dataset names.
        """
        self.dataset_lst = self.gi.show_matching_datasets(
            self.history_id, name_filter=name
        )
        self.ds_names = [k["name"] for k in self.dataset_lst]

    def get_histories(self):
        """Retrieves all histories.

        Returns:
            List: A list of histories.
        """
        self.histories = self.gi.histories.get_histories()
        return self.histories

    def clean_histories_for_display(self):
        """Cleans the histories for display.

        Returns:
            List: A list of cleaned histories.
        """
        self.histories_cleaned = [
            {
                "id": k["id"],
                "name": k["name"],
                "last_updated": k["update_time"],
            }
            for k in self.histories
            if k["deleted"] != True
        ]
        return self.histories_cleaned

    def set_tool(self, tool_id: str):
        """Sets the tool ID.

        Args:
            tool_id (str): The ID of the tool.
        """
        self.tool_id = tool_id
        self.logger.info(f"Tool Info: {self.gi.tools.show_tool(self.tool_id)}")

    def set_history(self, create: bool = True, hid: str = None, hname: str = None):
        """Sets the history.

        Args:
            create (bool, optional): Whether to create a new history. Defaults to True.
            hid (str, optional): The ID of the history. Defaults to None.
            hname (str, optional): The name of the history. Defaults to None.
        """
        # name is id or the name of the newly created history
        if create:
            if hname is None:
                hname = f"History created at {datetime.now()}"
            self.history_name = hname
            self.history_id = self.gi.histories.create_history(hname)["id"]
            self.logger.info(
                f"History Info: {self.gi.histories.show_history(self.history_id)}"
            )
        else:  # user wants to use the existing history
            self.history_id = hid
            self.history_name = self.gi.histories.show_history(self.history_id)["name"]

    def set_dataset(self, dataset_id: str):
        """Sets the dataset ID.

        Args:
            dataset_id (str): The ID of the dataset.
        """
        self.dataset_id = dataset_id
        self.logger.info(
            f"Dataset Info: {self.gi.datasets.show_dataset(self.dataset_id)}"
        )

    ########
    # Jobs #
    ########
    def show_job_status(self, job_id: str):
        """Shows the status of a job.

        Args:
            job_id (str): The ID of the job.
        """
        raise NotImplementedError

    ###########
    # Actions #
    ###########
    def upload_file(self, file_path: str):
        """Uploads a file.

        Args:
            file_path (str): The path to the file to upload.
        """
        upload = self.gi.tools.upload_file(file_path, self.history_id)
        self.set_dataset(upload["outputs"][0]["id"])
