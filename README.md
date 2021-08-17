# NuZTF
Python package for correlating ZTF data with external mutli-messenger triggers, created by [@robertdstein](https://github.com/robertdstein).
This package enables ZTF follow-up analysis of neutrinos/gravitational waves/gamma-ray bursts, built using the [AMPEL platform](https://arxiv.org/abs/1904.05922).

[![DOI](https://zenodo.org/badge/193068064.svg)](https://zenodo.org/badge/latestdoi/193068064) 
[![CI](https://github.com/robertdstein/ampel_followup_pipeline/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/robertdstein/ampel_followup_pipeline/actions/workflows/continous_integration.yml) 
![Binder](https://mybinder.org/badge_logo.svg) 
[![Coverage Status](https://coveralls.io/repos/github/desy-multimessenger/ampel_followup_pipeline/badge.svg?branch=master)](https://coveralls.io/github/desy-multimessenger/ampel_followup_pipeline?branch=master)

## Installation Instructions

The majority of required packages can be installed with the command:

```pip install -r requirements.txt```

You will need the IRSA login details with a ZTF-enabled account, to fully utilise all features.

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