from pathlib import Path

CONFIG_DIR = Path(__file__).parent.joinpath("config")

NUZTF_OUTPUT_DIR = Path.home().joinpath("Data/nuztf/")

SKYMAP_DIR = NUZTF_OUTPUT_DIR.joinpath(f"skymaps")
SKYMAP_DIR.mkdir(exist_ok=True, parents=True)
RESULTS_DIR = NUZTF_OUTPUT_DIR.joinpath("results")

CACHE_DIR = NUZTF_OUTPUT_DIR.joinpath(f"cache")

BASE_CANDIDATE_DIR = CACHE_DIR.joinpath("candidates")

CROSSMATCH_CACHE = CACHE_DIR.joinpath("crossmatch")
CROSSMATCH_CACHE.mkdir(exist_ok=True, parents=True)
CUTOUT_CACHE_DIR = CACHE_DIR.joinpath("cutouts")
CUTOUT_CACHE_DIR.mkdir(exist_ok=True, parents=True)
FLATPIX_CACHE_DIR = CACHE_DIR.joinpath("flatpix")
FLATPIX_CACHE_DIR.mkdir(exist_ok=True, parents=True)
