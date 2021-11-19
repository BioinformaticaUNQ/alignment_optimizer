import argparse
import logging
import sys

from project_alignment_optimizer import __version__
from project_alignment_optimizer.fib_model import fib

__author__ = "-"
__copyright__ = "-"
__license__ = "MIT"

_logger = logging.getLogger(__name__)


# controlador
def parse_args(args):
    parser = argparse.ArgumentParser(description="Just a Fibonacci demonstration")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="example_version: {version}".format(version=__version__),
    )
    parser.add_argument(
        dest="number",
        help="number-th Fibonacci number",
        type=int,
        metavar="INT"
    )
    parser.add_argument(
        "-vi",
        "--verbose_info",
        dest="loglevel",
        help="set loglevel to INFO",
        action="store_const",
        const=logging.INFO,
    )
    parser.add_argument(
        "-vd",
        "--verbose_debug",
        dest="loglevel",
        help="set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG,
    )
    return parser.parse_args(args)


def setup_logging(loglevel):
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S"
    )


# vista
def main(args):
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting crazy calculations...")
    print("The {}-th Fibonacci number is {}".format(args.number, fib(args.number)))
    _logger.info("Script ends here")

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run() 