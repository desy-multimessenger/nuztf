# AMPEL Follow-up Pipeline
Central repository for scripts handling ZTF Target-of-Opportunity requests, using the [AMPEL platform](https://arxiv.org/abs/1904.05922).

[![DOI](https://zenodo.org/badge/193068064.svg)](https://zenodo.org/badge/latestdoi/193068064) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/robertdstein/ZTF_Neutrino_ToO/master)[![CI](https://github.com/robertdstein/ampel_followup_pipeline/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/robertdstein/ampel_followup_pipeline/actions/workflows/continous_integration.yml)

## Installation Instructions

The majority of required packages can be installed with the command:

```pip install -r requirements.txt```

You will need the , as well as IRSA login details with a ZTF-enabled account, to fully utilise all features.

## Connecting to ampel

There are additional requirements to access ampel. To do so, you need to be a member of the ampel github group. Then you can additionally run:

```pip install -r ampel_requirements.txt```

In order to use the ampel tunnel, you need to have a separate shell that opens tunnels to the ampel database. You can do this using the following command in a terminal:

```ssh -L5432:localhost:5433 -L27020:localhost:27020 -L27018:localhost:27018 -L27026:localhost:27026 ztf-wgs.zeuthen.desy.de```

However, you need to make sure that your ~/.ssh/config file is configured directly for this to work.
