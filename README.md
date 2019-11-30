# ZTF_ToO_toolkit
Central repository for scripts handling ZTF neutrino Target of Opportunity requests

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/robertdstein/ZTF_Neutrino_ToO/master)

## Installation Instructions

The majority of required packages can be installed with the command:

```pip install -r requirements.txt```

You will need the , as well as IRSA login details with a ZTF-enabled account, to fully utilise all features.

## Connecting to ampel

There are additional requirements to access ampel. To do so, you need to be a member of the ampel github group. Then you can additionally run:

```pip install -r ampel_requirements.txt```

In order to use the ampel tunnel, you need to have a separate shell that opens a tunnel to the ampel database. You can do this using the following command in a terminal:

```ssh -L5432:localhost:5433 ztf-wgs.zeuthen.desy.de```

However, you need to make sure that your ~/.ssh/config file is configured directly for this to work.

In addition, you need a second port to connect to ampel resources and a third port to connect to the Catalog database 

THIS IS IN ADDITION TO THE PORT CONNECTING TO THE ARCHIVE DATABASE!

Open a new terminal, and into that terminal, run the following command:

```ssh -L27020:localhost:27020 ztf-wgs.zeuthen.desy.de```

THIS IS YET ANOTHER TUNNEL, NEEDED FOR USING EXTCAT TO QUERY CATALOGS

Open another new terminal (sorry), and within this terminal, run the following command:

```ssh -L27018:localhost:27018 ztf-wgs.zeuthen.desy.de```

