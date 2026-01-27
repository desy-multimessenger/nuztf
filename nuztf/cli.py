try:
    import typer
except ImportError:
    raise ImportError(
        "Please install typer if you want to use the CLI using `poetry install -E cli`"
    )

import logging
from typing import Annotated

from astropy.time import Time
from rich.console import Console
from rich.logging import RichHandler

from nuztf.neutrino_kafka_scanner import listen, scan_saved
from nuztf.neutrino_scanner import NeutrinoScanner

app = typer.Typer()
logger = logging.getLogger(__name__)


# --- Global callback (runs before every command) ---
@app.callback()
def main(
    ctx: typer.Context,
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        "-l",
        help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)",
        case_sensitive=False,
    ),
):
    """Global options for all nuztf commands."""
    # Normalize log level
    level = getattr(logging, log_level.upper(), None)
    if not isinstance(level, int):
        raise typer.BadParameter(f"Invalid log level: {log_level}")

    # Rich logging and printing
    logging.basicConfig(
        handlers=[RichHandler(rich_tracebacks=True, markup=True)],
    )
    logging.getLogger("nuztf").setLevel(level)
    console = Console()

    # Store log level in context for subcommands
    ctx.obj = {"log_level": level, "console": console}


@app.command()
def nu_classic(
    ctx: typer.Context,
    nu_name: Annotated[
        str, typer.Argument(..., help="Name of the neutrino, e.g. `IC200530A`")
    ],
    manual_info: Annotated[
        tuple[float, float, float, float, float, float, str],
        typer.Argument(
            ...,
            help="Provide manual position and time in the format "
            "RA, RA_LOWER, RA_UPPER, DEC, DEC_LOWER, DEC_UPPER, MJD",
        ),
    ] = None,
    gcn_filename: Annotated[
        str,
        typer.Option(
            "--gcn-filename",
            "-f",
            help="Filename to write GCN to, if None (default) print to console",
        ),
    ] = None,
):
    if manual_info is not None:
        manual_args = (
            nu_name,
            (manual_info[0], -manual_info[1], manual_info[2]),
            (manual_info[3], -manual_info[4], manual_info[5]),
            Time(manual_info[-1]),
        )
        nu = NeutrinoScanner(manual_args=manual_args)
    else:
        nu = NeutrinoScanner(nu_name=nu_name)
    nu.scan(console=ctx.obj["console"], gcn_filename=gcn_filename)


@app.command()
def nu_saved_kafka(
    ctx: typer.Context,
    nu_name: Annotated[
        str, typer.Argument(..., help="Name of the neutrino, e.g. `IC200530A`")
    ],
    gcn_filename: Annotated[
        str,
        typer.Option(
            "--gcn-filename",
            "-f",
            help="Filename to write GCN to, if None (default) print to console",
        ),
    ] = None,
):
    print("scanning")
    scan_saved(nu_name=nu_name, console=ctx.obj["console"], gcn_filename=gcn_filename)


@app.command()
def nu_listen(
    ctx: typer.Context,
    client_id: Annotated[str, typer.Argument(..., help="GCN Kafka Client ID")],
    client_secret: Annotated[str, typer.Argument(..., help="GCN Kafka Client Secret")],
    draft_directory: Annotated[
        str,
        typer.Option(
            "--draft-directory",
            "-d",
            help="Filename to write GCN to, if None (default) print to console",
        ),
    ] = None,
    replay: Annotated[
        str,
        typer.Option(
            "--from-utc-time",
            "-t",
            help="UTC time to start replaying from, in format YYYY-MM-DDTHH:MM:SS",
        ),
    ] = None,
):
    if replay is not None:
        try:
            replay = int(replay)
            logger.debug(f"Replaying last {replay} messages")
        except ValueError:
            logger.debug(f"Replaying from {replay}")

    listen(
        client_id=client_id,
        client_secret=client_secret,
        draft_directory=draft_directory,
        console=ctx.obj["console"],
        replay=replay,
    )


def cli_command():
    typer.run(main)
