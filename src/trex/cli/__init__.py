import logging
from pathlib import Path
import shutil

from trex.utils import NiceFormatter

logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    pass


def setup_logging(debug: bool) -> None:
    """
    Set up logging. If debug is True, then DEBUG level messages are printed.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(NiceFormatter())

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)


def add_file_logging(path: Path) -> None:
    file_handler = logging.FileHandler(path)
    root = logging.getLogger()
    root.addHandler(file_handler)


def make_output_dir(path, delete_if_exists):
    try:
        path.mkdir()
    except FileExistsError:
        if delete_if_exists:
            logger.debug(f'Re-creating folder "{path}"')
            shutil.rmtree(path)
            path.mkdir()
        else:
            raise
