import logging
import typer
from typing import Annotated
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

from nuztf.neutrino_scanner import NeutrinoScanner
from astropy.time import Time


def main(
        nu_name: str,
        logging_level:  Annotated[str, typer.Option("--logging-level", "-l")] = "INFO",
        gcn_filename: str | None = None
):
    logger = logging.getLogger("nuztf")
    logger.addHandler(RichHandler())
    logger.setLevel(logging_level)
    logger.propagate = False
    console = Console()
    console.print(f"Searching for candidates for {nu_name}", style="bold magenta")
    nu = NeutrinoScanner(nu_name)
    nu.query_ampel()
    nu.scan_area()
    nu.plot_overlap_with_observations(first_det_window_days=30.)
    jds = nu.observations.obsjd.unique()

    table = Table(title="Observations")
    table.add_column("Time", justify="right")
    table.add_column("Bands")
    table.add_column("Exp. Times")
    for jd in jds:
        m = nu.observations.obsjd == jd
        bands = nu.observations.band[m].unique()
        exp_times = nu.observations.exposure_time[m].unique()
        time = Time(jd, format="jd").ymdhms
        table.add_row(str(time), str(bands), str(exp_times))

    console.print(table)

    gcn_draft = nu.draft_gcn()
    if gcn_filename is not None:
        logger.info(f"Writing GCN to {gcn_filename}")
        with open(gcn_filename, "w") as f:
            f.write(gcn_draft)
    else:
        print(gcn_draft)


def cli_command():
    typer.run(main)
