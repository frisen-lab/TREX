"""
Simultaneous lineage TRacking and EXpression profiling of single cells using RNA-seq
"""
import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from . import cli as cli_package
from .cli import CommandLineError

from . import __version__

logger = logging.getLogger(__name__)


class HelpfulArgumentParser(ArgumentParser):
    """An ArgumentParser that prints full help on errors."""

    def __init__(self, *args, **kwargs):
        if "formatter_class" not in kwargs:
            kwargs["formatter_class"] = RawDescriptionHelpFormatter
        super().__init__(*args, **kwargs)

    def error(self, message):
        self.print_help(sys.stderr)
        args = {"prog": self.prog, "message": message}
        self.exit(2, "%(prog)s: error: %(message)s\n" % args)


def main(arguments=None):
    parser = HelpfulArgumentParser(description=__doc__, prog="trex")
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        subparser = subparsers.add_parser(
            module_name, help=module.__doc__.split("\n")[1], description=module.__doc__
        )
        subparser.set_defaults(func=module.main)
        module.add_arguments(subparser)

    args = parser.parse_args(arguments)
    subcommand = getattr(args, "func", None)
    if not subcommand:
        parser.error("Please provide the name of a subcommand to run")
    try:
        subcommand(args)
    except CommandLineError as e:
        logger.error(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
