# NuZTF
Python package for correlating ZTF data with external mutli-messenger triggers, created by [@robertdstein](https://github.com/robertdstein).
This package enables ZTF follow-up analysis of neutrinos/gravitational waves/gamma-ray bursts, built using the [AMPEL platform](https://arxiv.org/abs/1904.05922).

[![DOI](https://zenodo.org/badge/193068064.svg)](https://zenodo.org/badge/latestdoi/193068064)
[![CI](https://github.com/desy-multimessenger/nuztf/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/robertdstein/nuztf/actions/workflows/continous_integration.yml)
[![PyPI version](https://badge.fury.io/py/nuztf.svg)](https://badge.fury.io/py/nuztf)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/desy-multimessenger/nuztf/master)
[![Coverage Status](https://coveralls.io/repos/github/desy-multimessenger/nuztf/badge.svg?branch=master)](https://coveralls.io/github/desy-multimessenger/nuztf?branch=master)

## Installation Instructions

NuZTF can be directly installed with pip, giving the latest stable release:

```pip install nuztf```

Alternatively, the latest Github version of the code can be installed via pip:

```git clone https://github.com/desy-multimessenger/nuztf.git```

```pip install -e nuztf```

You will need the [IRSA login details](https://irsa.ipac.caltech.edu/account/signon/logout.do) with a ZTF-enabled account to fully utilise all features.

Additionally, you need an AMPEL API token. This can be obtained [here](https://ampel.zeuthen.desy.de/live/dashboard/tokens).

# Citing the code

If you make use of this code, please cite it! A DOI is provided by Zenodo, which can reference both the code repository, or specific releases:

[![DOI](https://zenodo.org/badge/193068064.svg)](https://zenodo.org/badge/latestdoi/193068064)

# Contributors

* Jannis Necker [@JannisNe](https://github.com/jannisne)
* Simeon Reusch [@simeonreusch](https://github.com/simeonreusch)
* Robert Stein [@robertdstein](https://github.com/robertdstein)

# Acknowledgements

This code stands on the shoulders of giants. We would particularly like to acknowledge:

* [Ampel](https://ampelproject.github.io/), created primarily by [@wombaugh](https://github.com/wombaugh), [@vbrinnel](https://github.com/vbrinnel) and [@jvansanten](https://github.com/jvansanten)
* [ztf_plan_obs](https://github.com/simeonreusch/ztf_plan_obs), created by [@simeonreusch](https://github.com/simeonreusch)
* [ztfquery](https://github.com/MickaelRigault/ztfquery), created by [@MickaelRigault](https://github.com/MickaelRigault)
