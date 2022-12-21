#!/usr/bin/env python3

try:
    import VarPredict.argP as argP
    import VarPredict.runVarPredict as runVarPredict
except ImportError:
    import argP
    import runVarPredict


def main():
    # set up and parse arguments
    parser = argP.buildParser()
    args = parser.parse_args()

    runVarPredict.run(args)

    return

if __name__ == "__main__":
    main()
