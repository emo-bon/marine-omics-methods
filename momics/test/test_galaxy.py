import pytest
import os

from unittest.mock import MagicMock, patch
from momics.galaxy.bio_blend import RemGalaxy


@pytest.fixture
def bbgalaxy():
    """
    Fixture that provides a Galaxy instance for testing.
    """
    with patch.dict(
        os.environ,
        {"GALAXY_URL": "http://fake-galaxy-instance", "GALAXY_API_KEY": "fake_api_key"},
    ):
        return RemGalaxy(url_var_name="GALAXY_URL", api_key_var_name="GALAXY_API_KEY")


def test_init(bbgalaxy):
    """
    Tests the initialization of the Galaxy instance.
    """
    assert (
        bbgalaxy.url == "http://fake-galaxy-instance"
    ), "The Galaxy URL should be set correctly"
    assert (
        bbgalaxy.api_key == "fake_api_key"
    ), "The Galaxy API key should be set correctly"


def test_get_histories(monkeypatch, bbgalaxy):
    """
    Tests the get_histories method.
    """
    mock_histories = [
        {"id": "1", "name": "History 1"},
        {"id": "2", "name": "History 2"},
    ]
    mock_get_histories = MagicMock(return_value=mock_histories)
    monkeypatch.setattr(bbgalaxy.gi.histories, "get_histories", mock_get_histories)

    histories = bbgalaxy.get_histories()

    assert histories == mock_histories, "The histories should be retrieved correctly"


# TODO: Too hard to fix now
# def test_create_history(monkeypatch, bcgalaxy):
#     """
#     Tests the set_history method for creating a new history.
#     """
#     mock_history = {"id": "3", "name": "New History"}
#     mock_create_history = MagicMock(return_value=mock_history)
#     monkeypatch.setattr(bcgalaxy.gi.histories, "create_history", mock_create_history)

#     bcgalaxy.set_history(create=True, hname="New History")

#     assert bcgalaxy.history_id == "3", "The new history ID should be set correctly"
#     assert (
#         bcgalaxy.history_name == "New History"
#     ), "The new history name should be set correctly"


def test_set_existing_history(monkeypatch, bbgalaxy):
    """
    Tests the set_history method for setting an existing history.
    """
    mock_history = {"id": "3", "name": "Existing History"}
    mock_show_history = MagicMock(return_value=mock_history)
    monkeypatch.setattr(bbgalaxy.gi.histories, "show_history", mock_show_history)

    bbgalaxy.set_history(create=False, hid="3")

    assert bbgalaxy.history_id == "3", "The existing history ID should be set correctly"
    assert (
        bbgalaxy.history_name == "Existing History"
    ), "The existing history name should be set correctly"


# def test_get_datasets_by_key(monkeypatch, bbgalaxy):
#     """
#     Tests the get_datasets_by_key method.
#     """
#     mock_datasets = [
#         {"name": "Dataset 1", "id": "value165464", "value": "val x"},
#         {"name": "Dataset 2", "id": "valueasdfasdf", "value": "val x"},
#     ]
#     mock_get_datasets = MagicMock(return_value=mock_datasets)
#     monkeypatch.setattr(bbgalaxy.gi.datasets, "get_datasets", mock_get_datasets)

#     datasets = bbgalaxy.get_datasets_by_key(key="id", value="value165464")

#     assert datasets == [("Dataset 1", "value165464")], "The datasets should be retrieved correctly by key and value"


# TODO: Too hard to fix now
# def test_upload_file(monkeypatch, bcgalaxy):
#     """
#     Tests the upload_file method.
#     """
#     mock_upload = {"outputs": [{"id": "dataset_id"}]}
#     mock_upload_file = MagicMock(return_value=mock_upload)
#     monkeypatch.setattr(bcgalaxy.gi.tools, "upload_file", mock_upload_file)

#     bcgalaxy.history_id = "history_id"
#     bcgalaxy.upload_file(file_path="fake_path")

#     assert (
#         bcgalaxy.dataset_id == "dataset_id"
#     ), "The dataset ID should be set correctly after upload"
