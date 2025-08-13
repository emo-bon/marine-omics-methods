import time
import json
import panel as pn

from momics.galaxy.bio_blend import RemGalaxy


class Gecco:
    """
    A class to interact with the GECCO tool in the Galaxy platform for comparative genomics.

    This class manages user authentication, history and dataset selection, file uploads,
    and submission of jobs to the GECCO tool via the Galaxy API.

    Args:
        params (dict): Dictionary of Panel widgets and parameters required for interaction.
    """

    def __init__(self, params):
        """
        Initialize the Gecco class with user parameters and Galaxy connection.

        Args:
            params (dict): Dictionary containing Panel widgets and configuration parameters.
        """
        self.exp = RemGalaxy("GALAXY_EARTH_URL", "GALAXY_EARTH_KEY")
        self.logger = self.exp.logger
        self.gecco_tool_id = "gecco"
        self.gecco_tool_name = "GECCO"
        self.gecco_tool_version = "0.1.0"
        self.gecco_tool_description = "GECCO is a tool for comparative genomics."
        self.select_history = params["select_history"]
        self.current_history_name = params["current_history_name"]
        self.current_history_id = params["current_history_id"]
        self.select_dataset = params["select_dataset"]
        self.current_file_name = params["current_file_name"]
        self.current_file_id = params["current_file_id"]
        self.file_source_checkbox = params["file_source_checkbox"]
        self.file_input = params["file_input"]
        self.mask = params["mask"]
        self.cds = params["cds"]
        self.threshold = params["threshold"]
        self.postproc = params["postproc"]
        self.antimash_sideload = params["antimash_sideload"]
        self.history_name = params["history_name"]
        self.logger.info("GECCO initialized.")
        self.debug = False

    def handle_login(self, clicks):
        """
        Handles user login and retrieves relevant data.

        Args:
            clicks (int): Number of times the login button has been clicked.
        """
        pn.state.notifications.info(f"User logged in: {self.exp.cfg.whoami()}")
        self.logger.info(f"You have clicked me {clicks} times")
        self.handle_get_histories(clicks)
        self.handle_get_datasets(clicks)

    def handle_get_histories(self, clicks):
        """
        Retrieves the user's Galaxy histories and updates the selection widget.

        Args:
            clicks (int): Number of times the login button has been clicked.
        """
        self.exp.get_histories()

        # clean histories dict for display, remove deleted histories and select fust relevant fields
        clean_histories = self.exp.clean_histories_for_display()
        self.select_history.options = clean_histories
        self.select_history.value = clean_histories[0]
        self.logger.info(f"{len(clean_histories)} histories found.")

        # update the current history name
        self.handle_update_current_history_name(self.select_history.value)
        self.handle_update_current_history_id(self.select_history.value)

    def handle_get_datasets(self, clicks):
        """
        Retrieves available datasets from Galaxy and updates the selection widget.

        Args:
            clicks (int): Number of times the login button has been clicked.
        """
        datasets, *_ = self.exp.get_datasets_by_key(
            "extension", ["fasta", "gbk", "embl", "genbank", "gb"]
        )
        self.logger.info(f"Datasets found: {datasets}")
        # fill the select_dataset widget
        self.select_dataset.options = datasets
        self.select_dataset.value = datasets[0]
        if self.debug:
            self.logger.info(
                f"Datasets selector options: {self.select_dataset.options}"
            )
            self.logger.info(
                f"Datasets selector value: {self.select_dataset.value}, {type(self.select_dataset.value)}"
            )
            self.logger.info(
                f"Datasets selector value[0]: {self.select_dataset.value[0]}"
            )
            self.logger.info(
                f"Datasets selector value[1]: {self.select_dataset.value[1]}"
            )
        self.handle_update_current_file_name(self.select_dataset.value)
        self.logger.info(f"{len(datasets)} datasets found.")

    def handle_update_current_file_name(self, value):
        """
        Updates the current file name and ID based on the selected dataset.

        Args:
            value (tuple): Tuple containing the file name and file ID.
        """
        self.current_file_name.value = value[0]
        self.current_file_id.value = value[1]

    def handle_update_current_history_name(self, value):
        """
        Updates the current history name based on the selected history.

        Args:
            value (dict): Dictionary containing history information.
        """
        self.current_history_name.value = value["name"]

    def handle_update_current_history_id(self, value):
        """
        Updates the current history ID based on the selected history.

        Args:
            value (dict): Dictionary containing history information.
        """
        self.current_history_id.value = value["id"]

    def handle_create_history(self, clicks):
        """
        Creates a new history in Galaxy if requested by the user.

        Args:
            clicks (int): Number of times the create history button has been clicked.
        """
        if clicks == 0:
            return
        if self.history_name.value != "":
            self.exp.set_history(hname=self.history_name.value)
        else:
            self.exp.set_history()

        # and update the select widget
        self.handle_get_histories(clicks)

    def handle_upload_dataset(self, clicks):
        """
        Uploads a dataset to Galaxy if the user chooses to upload from local source.

        Args:
            clicks (int): Number of times the upload button has been clicked.
        """
        if self.file_source_checkbox.value:
            pn.state.notifications.warning("You selected Galaxy source and not upload.")
        else:
            upload_data = self.exp.gi.tools.upload_file(
                self.file_input.value[0],
                self.select_history.value["id"],
            )

            # this needs to wait until the upload is done
            while (
                self.exp.gi.datasets.show_dataset(upload_data["outputs"][0]["id"])[
                    "state"
                ]
                != "ok"
            ):
                time.sleep(1)
            self.logger.info(f"Upload data: {upload_data}")
            uploaded_dataset_id = upload_data["outputs"][0]["id"]
            self.current_file_name.value = upload_data["outputs"][0]["name"]
            self.current_file_id.value = upload_data["outputs"][0]["id"]
            pn.state.notifications.success(f"Dataset {uploaded_dataset_id} uploaded.")
            self.logger.info(f"Dataset {uploaded_dataset_id} uploaded.")

    def handle_submit_gecco(self, clicks):
        """
        Submits a job to the GECCO tool in Galaxy with the selected parameters.

        Args:
            clicks (int): Number of times the submit button has been clicked.
        """
        if clicks == 0:
            pn.state.notifications.warning("You need to log in first.")
            return
        if self.current_file_id.value == "":
            pn.state.notifications.warning("No dataset selected.")
            self.logger.warning("No dataset selected.")
            return
        if self.current_history_id.value == "":
            pn.state.notifications.warning("No history selected.")
            self.logger.warning("No history selected.")
            return

        # create a submission dictionary
        submission_inputs = {
            "input": {
                "id": self.current_file_id.value,
                "src": "hda",
            },
            "mask": self.mask.value,
            "cds": self.cds.value,
            "threshold": self.threshold.value,
            "postproc": self.postproc.value,
            "antimash_sideload": self.antimash_sideload.value,
        }

        # Run the GECCO tool
        tool_run = self.exp.gi.tools.run_tool(
            history_id=self.current_history_id.value,
            tool_id=self.gecco_tool_id,
            tool_inputs=submission_inputs,
        )

        # Get job ID to monitor
        job_id = tool_run["jobs"][0]["id"]
        self.logger.info(f"GECCO tool job submitted with job ID: {job_id}")
        pn.state.notifications.info(f"GECCO tool job submitted with job ID: {job_id}")

        # add history ID and gecco ID to the submission_inputs
        submission_inputs["history_id"] = self.current_history_id.value
        submission_inputs["gecco_id"] = self.gecco_tool_id

        # save the job_id and submission details to local json file, name is the job_id
        with open(f"{job_id}.json", "w") as f:
            json.dump(submission_inputs, f)
        self.logger.info(f"Submission details saved to {job_id}.json")
