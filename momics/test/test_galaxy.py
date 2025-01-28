import pytest
import os

# import logging
# from datetime import datetime
from unittest.mock import MagicMock, patch
from momics.galaxy.blue_cloud import BCGalaxy


@pytest.fixture
def bcgalaxy():
    """
    Fixture that provides a BCGalaxy instance for testing.
    """
    with patch.dict(
        os.environ,
        {"GALAXY_URL": "http://fake-galaxy-instance", "GALAXY_API_KEY": "fake_api_key"},
    ):
        return BCGalaxy(url_var_name="GALAXY_URL", api_key_var_name="GALAXY_API_KEY")


def test_init(bcgalaxy):
    """
    Tests the initialization of the BCGalaxy instance.
    """
    assert (
        bcgalaxy.url == "http://fake-galaxy-instance"
    ), "The Galaxy URL should be set correctly"
    assert (
        bcgalaxy.api_key == "fake_api_key"
    ), "The Galaxy API key should be set correctly"


def test_get_histories(monkeypatch, bcgalaxy):
    """
    Tests the get_histories method.
    """
    mock_histories = [
        {"id": "1", "name": "History 1"},
        {"id": "2", "name": "History 2"},
    ]
    mock_get_histories = MagicMock(return_value=mock_histories)
    monkeypatch.setattr(bcgalaxy.gi.histories, "get_histories", mock_get_histories)

    histories = bcgalaxy.get_histories()

    assert histories == mock_histories, "The histories should be retrieved correctly"


def test_create_history(monkeypatch, bcgalaxy):
    """
    Tests the set_history method for creating a new history.
    """
    mock_history = {"id": "3", "name": "New History"}
    mock_create_history = MagicMock(return_value=mock_history)
    monkeypatch.setattr(bcgalaxy.gi.histories, "create_history", mock_create_history)

    bcgalaxy.set_history(create=True, hname="New History")

    assert bcgalaxy.history_id == "3", "The new history ID should be set correctly"
    assert (
        bcgalaxy.history_name == "New History"
    ), "The new history name should be set correctly"


def test_set_existing_history(monkeypatch, bcgalaxy):
    """
    Tests the set_history method for setting an existing history.
    """
    mock_history = {"id": "3", "name": "Existing History"}
    mock_show_history = MagicMock(return_value=mock_history)
    monkeypatch.setattr(bcgalaxy.gi.histories, "show_history", mock_show_history)

    bcgalaxy.set_history(create=False, hid="3")

    assert bcgalaxy.history_id == "3", "The existing history ID should be set correctly"
    assert (
        bcgalaxy.history_name == "Existing History"
    ), "The existing history name should be set correctly"


def test_get_datasets_by_key(monkeypatch, bcgalaxy):
    """
    Tests the get_datasets_by_key method.
    """
    mock_datasets = [
        {"name": "Dataset 1", "key": "value"},
        {"name": "Dataset 2", "key": "value"},
    ]
    mock_get_datasets = MagicMock(return_value=mock_datasets)
    monkeypatch.setattr(bcgalaxy.gi.datasets, "get_datasets", mock_get_datasets)

    datasets = bcgalaxy.get_datasets_by_key(key="key", value="value")

    assert datasets == [
        "Dataset 1",
        "Dataset 2",
    ], "The datasets should be retrieved correctly by key and value"


def test_upload_file(monkeypatch, bcgalaxy):
    """
    Tests the upload_file method.
    """
    mock_upload = {"outputs": [{"id": "dataset_id"}]}
    mock_upload_file = MagicMock(return_value=mock_upload)
    monkeypatch.setattr(bcgalaxy.gi.tools, "upload_file", mock_upload_file)

    bcgalaxy.history_id = "history_id"
    bcgalaxy.upload_file(file_path="fake_path")

    assert (
        bcgalaxy.dataset_id == "dataset_id"
    ), "The dataset ID should be set correctly after upload"
