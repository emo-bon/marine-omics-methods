import os
import panel as pn
from typing import List, Tuple
from pyngrok import ngrok
import pandas as pd
from IPython import get_ipython

from .utils import memory_load


def serve_app(template, env, name="panel app"):
    port = 4040
    while is_port_in_use(port):
        print("Port 4040 is in use, trying another port")
        port += 1
    print(f"Using port {port}")
    server = pn.serve(
        {name: template},
        port=port,
        address="127.0.0.1",
        threaded=True,
        websocket_origin="*",
    )

    if "google.colab" in str(get_ipython()) or env == "vscode":
        # server=pn.serve({"": template}, port=4040, address="127.0.0.1", threaded=True, websocket_origin="*")
        os.system(f"curl http://localhost:{port}")

        # Terminate open tunnels if exist
        ngrok.kill()

        # Setting the authtoken, get yours from https://dashboard.ngrok.com/auth
        NGROK_AUTH_TOKEN = os.getenv("NGROK_TOKEN")
        ngrok.set_auth_token(NGROK_AUTH_TOKEN)

        # Open an HTTPs tunnel on port 4040 for http://localhost:4040
        # if env == "vscode":
        public_url = ngrok.connect(addr=str(port))
        # else:
        #     public_url = ngrok.connect(port=str(port))

        print("Tracking URL:", public_url)
    else:
        pass
        # after development finished, this could be changed to np.serve()
        # server = pn.serve({name: template}, port=4040, address="127.0.0.1", threaded=True, websocket_origin="*")
        # template.servable()
    return server


def close_server(server, env):
    server.stop()
    if "google.colab" in str(get_ipython()) or env == "vscode":
        ngrok.disconnect(server)
        ngrok.kill()

    return


def is_port_in_use(port: int) -> bool:
    import socket

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("localhost", port)) == 0


def tax_finder_selector() -> Tuple[
    pn.widgets.Select,
    pn.widgets.Select,
    pn.widgets.TextInput,
    pn.widgets.Checkbox,
    pn.widgets.Checkbox,
]:
    select_table_tax = pn.widgets.Select(
        name="Taxonomic table",
        value="ssu",
        options=["ssu", "lsu"],
        description="Select a table for taxonomic search",
    )

    tax_level = pn.widgets.Select(
        name="Taxonomic level",
        value="all",
        options=[
            "all",
            "ncbi_tax_id",
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        description="Select a taxonomic search level",
    )

    search_term = pn.widgets.TextInput(
        name="Search term",
        value="",
        description="Enter a search term (string or NCBI tax ID)",
    )

    checkbox_exact_match = pn.widgets.Checkbox(
        name="Exact match of the search term",
        value=False,
    )

    log_scale_checkbox = pn.widgets.Checkbox(
        name="Log scale for abundance coloring",
        value=True,
    )

    return (
        select_table_tax,
        tax_level,
        search_term,
        checkbox_exact_match,
        log_scale_checkbox,
    )


def diversity_select_widgets(cat_columns: List[str], num_columns: List[str]) -> Tuple[
    pn.widgets.Select,
    pn.widgets.Select,
    pn.widgets.Select,
    pn.widgets.Select,
    pn.widgets.Select,
    pn.widgets.Checkbox,
]:
    """Creates selection widgets for alpha and beta diversity analysis.

    Args:
        cat_columns (List[str]): A list of categorical column names for alpha diversity.
        num_columns (List[str]): A list of numerical column names for beta diversity.

    Returns:
        Tuple[pn.widgets.Select, pn.widgets.Select, pn.widgets.Select, pn.widgets.Select, pn.widgets.Select]:
        A tuple containing selection widgets for alpha and beta diversity analysis.
    """
    select_table_alpha = pn.widgets.Select(
        name="Source table alphas",
        value="go",
        options=["go", "go_slim", "ips", "ko", "pfam"],
        description="Select a table for alpha diversity analysis",
    )

    select_cat_factor = pn.widgets.Select(
        name="Factor alpha",
        value=cat_columns[0],
        options=cat_columns,
        description="Categorical columns to compare alpha diversities",
    )

    select_table_beta = pn.widgets.Select(
        name="Source table beta",
        value="ssu",
        options=["ssu", "lsu"],
        description="Select a table for beta diversity analysis",
    )

    boptions = [
        "ncbi_tax_id",
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    select_taxon = pn.widgets.Select(
        name="Taxon",
        value=boptions[0],
        options=boptions,
        description="At which taxon level is beta diversity calculated",
    )

    full_columns = sorted(num_columns + cat_columns)
    select_factor_beta_all = pn.widgets.Select(
        name="Factor beta",
        value=full_columns[0],
        options=full_columns,
        description="Factor to visualize beta PCoA towards",
    )
    ## checkbox for beta matrix normalization
    checkbox_beta_norm = pn.widgets.Checkbox(
        name="Normalize beta matrix",
        value=False,
    )

    ret = (
        select_table_alpha,
        select_cat_factor,
        select_table_beta,
        select_taxon,
        select_factor_beta_all,
        checkbox_beta_norm,
    )
    return ret


#####################
# Create indicators #
#####################
def create_indicators_landing_page(
    df: pd.DataFrame,
) -> List[pn.indicators.Number]:
    """
    Generates a list of indicators for the landing page based on the provided aggregated DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame containing aggregated data with 'COMPLETED' and 'queued' columns.

    Returns:
        List[pn.indicators.Number]: A list of Panel indicators displaying the number of sequenced samples.
    """
    list_indicators = []

    nom1 = int(df["COMPLETED"]["filters"].sum())
    den1 = int(df["queued"]["filters"].sum()) + nom1
    list_indicators.append(
        pn.indicators.Number(
            name="Water samples sequenced",
            value=nom1,
            format="{value}" + f"/{den1}",
            width=150,
            font_size="34px",
            title_size="14px",
        ),
    )

    nom2 = int(df["COMPLETED"]["sediment"].sum())
    den2 = int(df["queued"]["sediment"].sum()) + nom2
    list_indicators.append(
        pn.indicators.Number(
            name="Sediment samples sequenced",
            value=nom2,
            format="{value}" + f"/{den2}",
            width=150,
            font_size="34px",
            title_size="14px",
        ),
    )

    try:
        nom3 = int(df["COMPLETED"]["water"].sum())
    except KeyError:
        nom3 = 0
    den3 = int(df["queued"]["filters_blank"].sum()) + nom3
    list_indicators.append(
        pn.indicators.Number(
            name="Blank Waters sequenced",
            value=nom3,
            format="{value}" + f"/{den3}",
            width=150,
            font_size="34px",
            title_size="14px",
        ),
    )

    try:
        nom4 = int(df["COMPLETED"]["sediment_blank"].sum())
    except KeyError:
        nom4 = 0
    den4 = int(df["queued"]["sediment_blank"].sum()) + nom4
    list_indicators.append(
        pn.indicators.Number(
            name="Blanks Sediments sequenced",
            value=nom4,
            format="{value}" + f"/{den4}",
            width=150,
            font_size="34px",
            title_size="14px",
        ),
    )
    return list_indicators


def create_indicators_diversity() -> (
    Tuple[pn.indicators.Progress, pn.indicators.Number]
):
    """Creates indicators for RAM usage.

    Returns:
        pn.FlexBox: A FlexBox containing RAM usage indicators.
    """
    used_gb, total_gb = memory_load()
    progress_bar = pn.indicators.Progress(
        name="Ram usage", value=int(used_gb / total_gb * 100), width=200
    )
    usage = pn.indicators.Number(
        value=used_gb,
        name="RAM usage [GB]",
        format="{value:,.2f}",
    )
    return progress_bar, usage
