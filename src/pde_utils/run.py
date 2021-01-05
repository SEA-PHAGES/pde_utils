import argparse
import timeit

from pde_utils.pipelines import (
                    build_pan, cluster_db, find_primers, make_db, pham_align)

VALID_PIPELINES = ["build_pan", "cluster_db", "find_primers",
                   "make_db", "pham_align"]


def main(unparsed_args):
    """Run a a pde_utils pipeline."""
    args = parse_args(unparsed_args)

    start = timeit.default_timer()

    if args.pipeline == "build_pan":
        build_pan.main(unparsed_args)
    elif args.pipeline == "cluster_db":
        cluster_db.main(unparsed_args)
    elif args.pipeline == "find_primers":
        find_primers.main(unparsed_args)
    elif args.pipeline == "make_db":
        make_db.main(unparsed_args)
    elif args.pipeline == "pham_align":
        pham_align.main(unparsed_args)

    stop = timeit.default_timer()

    print(f"\n\nPipeline completed.\nTime elapsed: {round(stop-start, 3)}s")


def parse_args(unparsed_args):
    """
    Use argparse to verify pipeline argument only.

    :param unparsed_args: raw_command line args
    :type unparsed_args: list
    :returns: ArgParse module parsed args.
    """

    RUN_HELP = """Command line script to call a pde_utils pipeline."""
    USAGE = """python3 -m pde_utils [pipeline]"""
    PIPELINE_HELP = """Name of the pde_utils pipeline to run."""

    parser = argparse.ArgumentParser(description=RUN_HELP, usage=USAGE)
    parser.add_argument("pipeline", type=str, choices=VALID_PIPELINES,
                        help=PIPELINE_HELP)

    args = parser.parse_args(unparsed_args[1:2])

    return args
