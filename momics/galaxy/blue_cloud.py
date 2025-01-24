import os
import logging
from datetime import datetime

from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.config import ConfigClient



class BCGalaxy:
    def __init__(self, url_var_name, api_key_var_name):
        self.logger = logging.getLogger('BCGalaxy')
        self.logger.setLevel(logging.DEBUG)
        # logging.debug("test")

        # url is a name of the variable saved in the .env file
        # the same for the api_key
        try:
            self.set_galaxy_env(url_var_name, api_key_var_name)
            self.cfg = ConfigClient(self.gi)
            self.logger.info(f'User: {self.cfg.whoami()}')
            self.logger.info(f'Galaxy version: {self.cfg.get_version()}')
        except TypeError:
            self.logger.error("Please provide valid Galaxy URL and API key!")
        except Exception as e:
            self.logger.error(f"Error: {e}")
    
    ###########
    # Getters #
    ###########
    def set_galaxy_env(self, url_var_name, api_key_var_name):
        self.url = os.getenv(url_var_name)
        self.api_key = os.getenv(api_key_var_name)
        self.gi = GalaxyInstance(self.url, key=self.api_key)

    def get_datasets_by_key(self, key, value):
        lst_dict = [k for k in self.gi.datasets.get_datasets() if 
                    key in k and
                    k[key] == value]

        self.ds_names = [k["name"] for k in lst_dict]
        return self.ds_names
    
    # def get_datasets_id(self, name):
    #     files = [k["id"] for k in self.gi.datasets.get_datasets() if k["name"] == name]
    #     if len(files) > 1:
    #         self.logger.warning(f"Multiple datasets with the same name: {name}")
    #     self.dataset_id = files[0]
    #     return self.dataset_id

    def get_datasets(self, name=None):
        self.dataset_lst = self.gi.show_matching_datasets(self.history_id, name_filter=name)
        self.ds_names = [k["name"] for k in self.dataset_lst]
    
    def get_histories(self):
        self.histories = self.gi.histories.get_histories()

    ###########
    # Setters #
    ###########
    def set_tool(self, tool_id):
        self.tool_id = tool_id
        self.logger.info(f"Tool Info: {self.gi.tools.show_tool(self.tool_id)}")

    def set_history(self, create: bool = True, hid=None, hname=None):
        # name is id or the name of the newly created history
        if create:
            if hname is None:
                hname = f"History created at {datetime.now()}"
            self.history_name = hname
            self.history_id = self.gi.histories.create_history(hname)['id']
            self.logger.info(
                f"History Info: {self.gi.histories.show_history(self.history_id)}")
        else:  # user wants to use the existing history
            self.history_id = hid
            self.history_name = self.gi.histories.show_history(self.history_id)['name']

    def set_dataset(self, dataset_id):
        self.dataset_id = dataset_id
        self.logger.info(f"Dataset Info: {self.gi.datasets.show_dataset(self.dataset_id)}")

    ########
    # Jobs #
    ########
    def show_job_status(self, job_id):
        raise NotImplementedError

    ###########
    # Actions #
    ###########
    def upload_file(self, file_path):
        upload = self.gi.tools.upload_file(file_path, self.history_id)
        self.set_dataset(upload['outputs'][0]['id'])