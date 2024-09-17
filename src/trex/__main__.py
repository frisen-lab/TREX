"""
Simultaneous lineage TRacking and EXpression profiling of single cells using RNA-seq
"""

import ast
import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from . import cli as cli_package
from .cli import CommandLineError, setup_logging

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
    subcommand_name = get_subcommand_name(arguments)
    module = importlib.import_module("." + subcommand_name, cli_package.__name__)
    parser = HelpfulArgumentParser(description=__doc__, prog="trex")
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="Print some extra debugging messages",
    )

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser(
        subcommand_name,
        help=module.__doc__.strip().split("\n", maxsplit=1)[0],
        description=module.__doc__,
    )
    module.add_arguments(subparser)
    args = parser.parse_args(arguments)
    setup_logging(args.debug)

    del args.debug
    try:
        module.main(args)
    except CommandLineError as e:
        logger.error("trex error: %s", str(e))
        logger.debug("Command line error. Traceback:", exc_info=True)
        sys.exit(1)


def get_subcommand_name(arguments) -> str:
    """
    Parse arguments to find out which subcommand was requested.

    This sets up a minimal ArgumentParser with the correct help strings.

    Because help is obtained from a module’s docstring, but importing each module
    makes startup slow, the modules are only parsed with the ast module and
    not fully imported at this stage.

    Return:
        subcommand name
    """
    parser = HelpfulArgumentParser(description=__doc__, prog="trex")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    for module_name, docstring in cli_modules(cli_package):
        help = docstring.strip().split("\n", maxsplit=1)[0].replace("%", "%%")
        subparser = subparsers.add_parser(
            module_name, help=help, description=docstring, add_help=False
        )
        subparser.set_defaults(module_name=module_name)
    args, _ = parser.parse_known_args(arguments)
    module_name = getattr(args, "module_name", None)
    if module_name is None:
        parser.error("Please provide the name of a subcommand to run")
    return module_name


def cli_modules(package):
    """
    Yield (module_name, docstring) tuples for all modules in the given package.
    """
    modules = pkgutil.iter_modules(package.__path__)
    for module in modules:
        spec = importlib.util.find_spec(package.__name__ + "." + module.name)
        with open(spec.origin) as f:
            mod_ast = ast.parse(f.read())
        docstring = ast.get_docstring(mod_ast, clean=False)
        yield module.name, docstring


if __name__ == "__main__":
    main()
