#!/usr/bin/env python3

from __future__ import print_function

import argparse
import os
import sys

from time import gmtime, strftime
from ruamel.yaml import YAML
try:
    from termcolor import colored
except ImportError:
    sys.stderr.write("Warning: The termcolor package is unavailable\n")
    sys.stderr.write("         The --fancy option will be ignored\n")
    def colored(text, fgcolor, bgcolor=None, **kwargs):
        return text


                    # ------------------------------------ #


# Color attributes for colorized reports
fancy_color = {
    "????": ["grey", "on_yellow", "bold,dark"],
    "enfl": ["grey", "on_yellow", "bold,dark"],
    "eout": ["yellow", "on_red", "bold"],
    "eref": ["white", "on_red", "bold"],
    "etol": ["yellow", None, "bold"],
    "xtfl": ["yellow", "on_green", "bold"],
    "pass": ["green", None, "bold"],
    "skip": ["grey", None, "bold"],
}

# Test status formatter for both BW and color reports
def fancy_status(stat_str, colorize=False):

    # Maximum status length is 4
    if ( len(stat_str) >= 4 ):
        fancy_str = stat_str[0:4]
    else:
        fancy_str = stat_str

    # Minimum status length is 4
    if ( len(fancy_str) == 0 ):
        fancy_str = "????"
    if ( len(fancy_str) <= 3 ):
        fancy_str += " "
    if ( len(fancy_str) <= 3 ):
        fancy_str = " " + fancy_str
    if ( len(fancy_str) <= 3 ):
        fancy_str += " "

    # Colorize output when requested
    if ( colorize ):
        if ( fancy_str in fancy_color ):
            fgcolor, bgcolor, deco = fancy_color[fancy_str]
            if ( len(deco) > 0 ):
                deco = deco.split(",")
            else:
                deco = []
        else:
            fgcolor = "yellow"
            bgcolor = "on_blue",
            deco = ["bold"]
        if ( bgcolor ):
            fancy_str = colored(fancy_str.upper(), fgcolor, bgcolor, attrs=deco)
        else:
            fancy_str = colored(fancy_str.upper(), fgcolor, attrs=["bold"])
        fancy_str = colored("[", "white") + fancy_str + colored("]", "white")
    else:
        fancy_str = "[" + fancy_str.upper() + "]"

    return fancy_str


                    # ------------------------------------ #


# Process command-line arguments
parser = argparse.ArgumentParser(
    description="Checks SIESTA test results against references")
parser.add_argument("-a", "--abs-tol", type=float, default=1.0e-7, 
    help="Absolute tolerance for comparisons (overrides configuration)")
parser.add_argument("-c", "--config", default="siesta-testsuite.yml",
    help="Config file containing the SIESTA test suite specifications")
parser.add_argument("-d", "--details", action="store_true", default=False,
    help="Display detailed information about each test")
parser.add_argument("-f", "--fancy", action="store_true", default=False,
    help="Print fancy colorized report (for terminals supporting it)")
parser.add_argument("-o", "--output", default="tests-report.yml",
    help="File to store a YAML version of the report with full details")
parser.add_argument("-r", "--refdir", default="YAML_Refs",
    help="Directory containing reference YAML files")
parser.add_argument("-t", "--testdir", default=".",
    help="Directory containing the test cases with their outputs")
parser.add_argument("-y", "--yaml", default=os.path.join("work", "OUTVARS.yml"),
    help="Relative path where each test stores its YAML data")
args = parser.parse_args()

# Check command-line arguments
if ( not os.path.exists(args.config) ):
    parser.error("config file not found: '%s'" % args.config)
if ( not os.path.isdir(args.refdir) ):
    parser.error("reference directory not found: '%s'" % args.refdir)
if ( not os.path.isdir(args.testdir) ):
    parser.error("tests directory not found: '%s'" % args.testdir)

# Banner
print("""\
========================================
Summary of the SIESTA Test Suite results
========================================

Introduction
------------

Each relevant test case of the SIESTA Test Suite has been run and its results
compared to a reference file.

In the following, the result of each comparison is represented in a compact
way at the beginning of each line through a 4-letter keyword surrounded by
square brackets. It is followed by a description of the correpsonding test,
condensed by default, verbose if requested.

Here is the meaning of each keyword:

- ENFL: The test case succeeded while it was expected to fail.
- EOUT: Error reading YAML output file of the test case.
- EREF: Error reading YAML reference file of the test case.
- ETOL: Some values are beyond permitted tolerances.
- PASS: The test case meets all requirements.
- SKIP: The test case has not been executed.
- XTFL: The test failed as expected.

A test is successful only if its result is "PASS".

A report with full details of the tests is stored in the output file listed
below in YAML format.


Test suite parameters
---------------------

The test comparison script has been run with the following parameters:

- Workdir       : %s
- Configuration : %s
- Reference dir : %s
- Test dir      : %s
- Output file   : %s


Results
-------
""" % (os.getcwd(), args.config, args.refdir, args.testdir, args.output))

# Load global configuration
yaml_doc = YAML()
yaml_cfg = {"tests": []}
with open(args.config, "r") as cfg_file:
    yaml_cfg = yaml_doc.load(cfg_file)

# Perform test comparisons
siesta_tests = []
siesta_index = {}
siesta_xfail = []
for tcase in yaml_cfg["tests"]:

    # Init
    siesta_index[tcase["name"]] = len(siesta_tests)
    siesta_tests.append({"name": tcase["name"], "title": tcase["title"]})
    if ( ("fail_expect" in tcase) and (tcase["fail_expect"] == "yes") ):
        siesta_xfail.append(tcase["name"])

    # Check that test reference exists
    ref_path = os.path.join(args.refdir, "%s.yml" % tcase["name"])
    if ( not os.path.exists(ref_path) ):
        siesta_tests[-1]["result"] = "eref"
        siesta_tests[-1]["message"] = "missing reference file"
        continue

    # Load test reference
    yaml_doc = YAML()
    with open(ref_path, "r") as ref_file:
        yaml_ref = yaml_doc.load(ref_file)

    # Check that test output exists
    out_path = os.path.join(args.testdir, tcase["name"], args.yaml)
    if ( not os.path.exists(out_path) ):
        siesta_tests[-1]["result"] = "eout"
        siesta_tests[-1]["message"] = "missing output file"
        continue

    # Load test output
    yaml_doc = YAML()
    with open(out_path, "r") as out_file:
        yaml_out = yaml_doc.load(out_file)

    # Compare energies
    siesta_tests[-1]["tolerances"] = {}
    ref_vars = yaml_ref["energies"]
    out_vars = yaml_out["energies"]
    tc_good = []
    tc_fail = []
    tc_skip = []
    for (key, val) in ref_vars.items():
        if ( ("tolerances" in tcase.keys()) and (key in tcase["tolerances"]) ):
            tc_tol = tcase["tolerances"][key]
        elif ( key in yaml_cfg["tolerances"] ):
            tc_tol = yaml_cfg["tolerances"][key]
        else:
            tc_tol = args.abs_tol
        siesta_tests[-1]["tolerances"][key] = tc_tol

        if ( (key in out_vars) and (abs(out_vars[key] - val) < tc_tol) ):
            tc_good.append(key)
        else:
            tc_fail.append(
              {"name": key, "value": out_vars[key], "expected": val})

    # Store comparison results
    tc_skip = sorted([item for item in out_vars.keys() if not item in ref_vars.keys()])
    siesta_tests[-1]["good"] = tc_good
    siesta_tests[-1]["fail"] = tc_fail
    siesta_tests[-1]["skip"] = tc_skip

    # Store final test result
    if ( len(tc_fail) == 0 ):
        siesta_tests[-1]["result"] = "pass"
        siesta_tests[-1]["message"] = "all values within tolerance"
    else:
        siesta_tests[-1]["result"] = "etol"
        siesta_tests[-1]["message"] = "some values beyond tolerance"

# Write down full report
yaml_doc = YAML()
yaml_doc.default_flow_style=False
yaml_doc.indent(mapping=4, sequence=4, offset=2)
report = {"siesta_tests": siesta_tests,
          "defaults": {"tolerance": args.abs_tol}}
with open(args.output, "w") as rep_file:
    rep_file.write("%YAML 1.2\n---\n\n")
    yaml_doc.dump(report, rep_file)
    rep_file.write("\n...\n")

# Display test results
for tcase in siesta_tests:

    if ( tcase["name"] in siesta_xfail ):
        if ( tcase["result"] == "pass" ):
            tcase["result"] = "enfl"
        else:
            tcase["result"] = "xtfl"
        tcase["message"] += " (expected to fail)"

    if ( not args.details ):
        titlen = 65 - len(tcase["name"])
        if ( len(tcase["title"]) > titlen ):
            tcase["title"] = tcase["title"][:titlen] + "..."

    if ( args.fancy ):
        tcase["name"] = colored(tcase["name"], "white", attrs=["bold"])
        tcase["title"] = colored(tcase["title"], "white", attrs=["dark"])

    print("-", fancy_status(tcase["result"], args.fancy), tcase["name"],
        "> %s" % tcase["title"])

    if ( args.details ):
        if ( not tcase["result"] in ["eout", "eref"] ):
            print("")
            print("  - Energies/good: %s", tcase["good"])
            print("  - Energies/fail: %s", tcase["fail"])
            print("  - Energies/skip: %s", tcase["skip"])
            print("")

# Display statistics
ntests = len(siesta_tests)
ntpass = len([tcase for tcase in siesta_tests if tcase["result"] == "pass"])
print("""\


Statistics
----------

Global results of the test suite:

- number of tests  : %4d
- successful tests : %4d (%5.1f%%)
- failed tests     : %4d (%5.1f%%)

""" % (ntests, ntpass, ntpass*100.0/ntests, ntests-ntpass, (ntests-ntpass)*100.0/ntests))

# Display footer
print("""\
.. note::

   This report is valid ReStructuredText, unless you ask for colorized output.
   You can render it in HTML, PDF and other formats using the Docutils package.
   Please consult `http://docutils.sourceforge.net/`_ for details.

Document generated on %s (UTC).
""" % strftime("%Y/%m/%d %H:%M:%S +0000",gmtime()))
