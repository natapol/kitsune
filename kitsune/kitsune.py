#!/usr/bin/env python3
"""
KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics
"""

import importlib
import sys

__author__ = "Natapol Pornputtapong (natapol.p@chula.ac.th)"
__version__ = "1.3.5"
__date__ = "May 16, 2024"

# Define the set of commands with help messages
COMMANDS = {
    "acf": "Compute average number of common features between signatures",
    "cre": "Compute cumulative relative entropy",
    "dmatrix": "Compute distance matrix",
    "kopt": ("Compute recommended choice (optimal) of kmer within a given "
             "kmer interval for a set of genomes using the cre, acf and ofc"),
    "ofc": "Compute observed feature frequencies"
}


def help():
    """
    Print the list of available commands and exit
    """

    print("\nusage: kitsune <command> [<args>]\n")

    print("Available commands:")

    for cmd in COMMANDS:
        print("\t{}\t{}".format(cmd, COMMANDS[cmd]))
    
    print(
        ("\nUse --help in conjunction with one of the commands above for a list of available options "
         "(e.g. kitsune {} --help)").format(list(COMMANDS.keys())[0])
    )


def main():
    if len(sys.argv) < 2:
        # Print --help in case of no arguments
        help()
        sys.exit()

    # Get command
    cmd = sys.argv[1]

    if cmd == "--version" or cmd == "-v":
        print("kitsune v{} ({})".format(__version__, __date__))
        sys.exit()
    
    if cmd == "--help" or cmd == "-h":
        help()
        sys.exit()

    # Check whether cmd is a valid command
    if cmd not in COMMANDS:
        raise Exception(
            'Command "{}" not found! Type "kitsune --help" for a list of available commands'.format(sys.argv[1])
        )

    # Remove "kitsune" and command from sys.argv
    sys.argv.remove(sys.argv[1])
    sys.argv.remove(sys.argv[0])

    # Import command and run
    module = importlib.import_module("kitsune.modules.{}".format(cmd))
    module.run(sys.argv)


if __name__ == "__main__":
    main()
