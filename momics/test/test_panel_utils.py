import pytest
import panel as pn
from momics.panel_utils import diversity_select_widgets


@pytest.fixture
def sample_columns():
    """
    Fixture that provides sample categorical and numerical columns for testing.
    """
    cat_columns = ["factor1", "factor2", "factor3"]
    num_columns = ["num1", "num2", "num3"]
    return cat_columns, num_columns


def test_diversity_select_widgets(sample_columns):
    """
    Tests the diversity_select_widgets function.
    """
    cat_columns, num_columns = sample_columns
    widgets = diversity_select_widgets(cat_columns, num_columns)

    # Check if the result is a tuple
    assert isinstance(widgets, tuple), "The result should be a tuple"

    # Check if the tuple contains 5 elements
    assert len(widgets) == 6, "The tuple should contain 5 elements"

    # Check if each element in the tuple is a pn.widgets.Select instance
    for widget in widgets[:-1]:
        assert isinstance(
            widget, pn.widgets.Select
        ), f"Expected pn.widgets.Select, but got {type(widget)}"

    # Check if the last element in the tuple is a pn.widgets.Checkbox instance
    assert isinstance(
        widgets[-1], pn.widgets.Checkbox
    ), f"Expected pn.widgets.Checkbox, but got {type(widgets[-1])}"

    # Check if the widgets have the correct initial values and options
    assert (
        widgets[0].value == "go"
    ), "The initial value of the first widget should be 'go'"
    assert widgets[0].options == [
        "go",
        "go_slim",
        "ips",
        "ko",
        "pfam",
    ], "The options of the first widget are incorrect"

    assert (
        widgets[1].value == cat_columns[0]
    ), f"The initial value of the second widget should be '{cat_columns[0]}'"
    assert (
        widgets[1].options == cat_columns
    ), "The options of the second widget are incorrect"

    assert (
        widgets[2].value == "ssu"
    ), "The initial value of the third widget should be 'ssu'"
    assert widgets[2].options == [
        "ssu",
        "lsu",
    ], "The options of the third widget are incorrect"

    assert (
        widgets[3].value == "ncbi_tax_id"
    ), "The initial value of the fourth widget should be 'ncbi_tax_id'"
    assert widgets[3].options == [
        "ncbi_tax_id",
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ], "The options of the fourth widget are incorrect"

    full_columns = sorted(num_columns + cat_columns)
    assert (
        widgets[4].value == full_columns[0]
    ), f"The initial value of the fifth widget should be '{full_columns[0]}'"
    assert (
        widgets[4].options == full_columns
    ), "The options of the fifth widget are incorrect"
