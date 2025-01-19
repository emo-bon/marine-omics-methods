import panel as pn
from .helpers import memory_load


def create_widgets(cat_columns, num_columns, styles):
    select_table_alpha = pn.widgets.Select(
        name="Source table alphas",
        value="go",
        options=["go", "go_slim", "ips", "ko", "pfam"],
        description="Select a table for alpha diversity analysis",
    )

    # TODO: probably not all of these categories make sense
    select_cat_factor = pn.widgets.Select(
        name="Factor alpha",
        value=cat_columns[0],
        options=cat_columns,
        description="Categorical columns to compare alpha diversities",
    )

    # For Beta diversity
    select_table_beta = pn.widgets.Select(
        name="Source table beta",
        value="SSU",
        options=["SSU", "LSU"],
        description="Select a table for beta diversity analysis",
    )

    boptions = ["ncbi_tax_id", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]  # Options for table selection
    select_taxon = pn.widgets.Select(
        name="Taxon",
        value=boptions[0],
        options=boptions,
        description="At which taxon level is beta diversity calculated",
    )

    select_factor_beta = pn.widgets.Select(
        name="Factor beta",
        value=num_columns[0],
        options=num_columns,
        description="Factor to visualize beta PCoA towards",
    )

    # indicators
    # used_gb, total_gb = memory_load()
    # indicators = pn.FlexBox(
    #     pn.indicators.Progress(
    #         name='Ram usage',
    #         value=int(used_gb / total_gb * 100),
    #         width=200
    #     ),
    #     pn.indicators.Number(
    #         value=used_gb, name="RAM usage [GB]",
    #         format="{value:,.1f}",
    #         styles=styles
    #     ),
    #     pn.indicators.Number(
    #         value=10, name="Not implemented", format="{value:,.0f}", styles=styles
    #     ),
    # )
    ret = (
        select_table_alpha,
        select_cat_factor,
        select_table_beta,
        select_taxon,
        select_factor_beta,
        # indicators,
    )
    return ret


def create_indicators():
    used_gb, total_gb = memory_load()
    indicators = pn.FlexBox(
        pn.indicators.Progress(
            name='Ram usage',
            value=int(used_gb / total_gb * 100),
            width=200
        ),
        pn.indicators.Number(
            value=used_gb, name="RAM usage [GB]",
            format="{value:,.1f}",
        ),
    )
    return indicators