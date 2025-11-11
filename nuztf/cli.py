try:
    import typer
except ImportError:
    raise ImportError(
        "Please install typer if you want to use the CLI using `poetry install -E cli`"
    )

import logging
from pathlib import Path
from typing import Annotated

from rich.console import Console
from rich.logging import RichHandler

from nuztf.neutrino_scanner import NeutrinoScanner


def main(
    nu_name: Annotated[
        str, typer.Argument(..., help="Name of the neutrino, e.g. `IC200530A`")
    ],
    logging_level: Annotated[str, typer.Option("--logging-level", "-l")] = "INFO",
    gcn_filename: Annotated[
        str | Path,
        typer.Option(
            "--gcn-filename",
            "-f",
            help="Filename to write GCN to, if None (default) print to console",
        ),
    ] = None,
    rich_handler: Annotated[
        bool, typer.Option("--rich", "-r", help="Use a nice logging markup")
    ] = True,
    stream_handler: Annotated[
        bool,
        typer.Option("--stream", "-s", help="Use the standard stdout logging handler"),
    ] = False,
):
    logger = logging.getLogger("nuztf")
    if rich_handler:
        logger.addHandler(RichHandler())
    if not stream_handler:
        logger.propagate = False
    logger.setLevel(logging_level)
    console = Console()
    console.print(f"Searching for candidates for {nu_name}", style="bold magenta")
    nu = NeutrinoScanner(nu_name)
    nu.scan(console=console, gcn_filename=gcn_filename)


def cli_command():
    typer.run(main)
