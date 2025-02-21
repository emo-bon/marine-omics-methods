import os
import panel as pn
from typing import List, Tuple
from pyngrok import ngrok
from IPython import get_ipython

from .utils import memory_load


def serve_app(template, env, name="panel app"):
    port = 4040
    while is_port_in_use(port):
        print("Port 4040 is in use, trying another port")
        port += 1
    print(f"Using port {port}")
    server=pn.serve({name: template}, port=port, address="127.0.0.1", threaded=True,
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
        if env == "vscode":
            public_url = ngrok.connect(addr=str(port))
        else:
            public_url = ngrok.connect(port=str(port))
        
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
        return s.connect_ex(('localhost', port)) == 0


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
        value="SSU",
        options=["SSU", "LSU"],
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


def create_indicators() -> Tuple[pn.indicators.Progress, pn.indicators.Number]:
    """Creates indicators for RAM usage.

    Returns:
        pn.FlexBox: A FlexBox containing RAM usage indicators.
    """
    used_gb, total_gb = memory_load()
    # indicators = pn.FlexBox(
    progress_bar = pn.indicators.Progress(
            name="Ram usage", value=int(used_gb / total_gb * 100), width=200
        )
    usage = pn.indicators.Number(
            value=used_gb,
            name="RAM usage [GB]",
            format="{value:,.2f}",
        )
    return progress_bar, usage
