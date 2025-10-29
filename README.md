# marine-omics

[![Python Version](https://img.shields.io/pypi/pyversions/momics.svg?color=green)](https://python.org)
![PyPI - Version](https://img.shields.io/pypi/v/marine-omics)
[![Read the Docs](https://img.shields.io/readthedocs/marine-omics)](https://marine-omics-methods.readthedocs.io/en/latest/)
[![tests](https://github.com/palec87/marine-omics/workflows/tests/badge.svg)](https://github.com/palec87/marine-omics/actions)
[![codecov](https://codecov.io/gh/emo-bon/marine-omics-methods/branch/main/graph/badge.svg)](https://codecov.io/gh/emo-bon/marine-omics-methods)
[![PyPI Downloads](https://static.pepy.tech/badge/marine-omics/month)](https://pepy.tech/projects/marine-omics)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Package of utilities for FAIR-Ease demo workflows.

These methods are currently use only for FAIR-EASE [pilot demos](https://github.com/palec87/momics-demos), but eventually they can serve for general purpose manipulation of metagenomic data, locally and in VREs.

The idea is to provide testable methods to allow as much flexibility and remixing of the functionalities provided.

## Installation

One of the dependencies is currently not released to PyPI and is in active development, therefore you need to manually install it before `marine-omics` itself:

```bash
# UDAL data query layer
pip install git+https://github.com/fair-ease/py-udal-mgo.git

pip install marine-omics
```

## European Marine Omics Biodiversity Observation Network

Specifically, we aim primarily to manipulate EMO-BON marine genomics sampling data and metadata from ENA project [PRJEB51688](https://www.ebi.ac.uk/ena/browser/view/PRJEB51688). The interactive dashboards and jupyter notebooks built on top of this repository can be found [here](https://github.com/emo-bon/momics-demos/tree/main).

The methods are mixture of statistical methods, plotting functionalities, metadata and data handling utilities and generators of `holoviz panel` widgets and panes. Experimental integration to `Galaxy` uses a wrapper around [bioblend](https://bioblend.readthedocs.io/en/latest/).
