
# Copyright 2014 Fabio Cascioli, Jonas Lindert, Philipp Maierhoefer, Stefano Pozzorini
#
# This file is part of OpenLoops.
#
# OpenLoops is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenLoops is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenLoops.  If not, see <http://www.gnu.org/licenses/>.


import sys

if sys.version_info[:2] < (2,7):
    print "Python 2.7 required."
    sys.exit(1)

import os
import openloops_legacy as openloops
import keywordoptions

false_values = ("false", "f", "none", "no", "n", "off", "0")
true_values = ("true", "t", "all", "yes", "y", "on", "1")

args = sys.argv
check = False

if "--check" in args:
    args = [arg for arg in args if arg != "--check"]
    check = True

ko = keywordoptions.KeywordOptions(strict = False)

ko.add("energy", converter=float)
ko.add("subprocess", converter=openloops.parse_subprocess_argument)
ko.add("n", converter=int)
ko.add("verbose", converter=int, default=1)
#ko.add(openloops.Parameters.names, group="parameters")
ko.add("save", default="False")
if check:
    ko.add("data", default="True")


def strip_stringlist(ls):
    """
    Take a list of strings, cut away comments starting with "#",
    strip off whitespace, and remove empty strings.
    """
    strippedls = [li.split("#")[0].strip() for li in ls]
    strippedls = [li for li in strippedls if li != ""]
    return strippedls

def split_first(s):
    ls = s.split(" ", 1)
    if len(ls) == 1:
        ls.append("")
    return ls


def default_data_file(process_name):
    return os.path.join(openloops.test_data_folder, process_name + ".ptd")


def data_files(data, save, process_name=None):
    # data is None/false/0 -> no src, or true/1 -> auto src, or a file name
    if not data or data.lower() in false_values:
        src = None
    elif data.lower() in true_values:
        src = default_data_file(process_name)
    else:
        src = data
    # save is false/0/new -> no dst, or true/1 -> auto dst = src if defined, otherwise default dst.
    # save=new will save to src if it doesn't exist.
    if save.lower() in false_values or save.lower() == "new":
        dst = None
    elif save.lower() in true_values:
        if src is None:
            dst = default_data_file(process_name)
        else:
            dst = src
    else:
        dst = save
    return src, dst


try:
    with open(args[1]) as fh:
        file_input = True
        test_definitions = [split_first(td) for td in strip_stringlist(fh.readlines())]
except IOError:
    file_input = False
    test_definitions = [[args[1], args[2:]]]


if file_input:
    global_options = [td[1] for td in test_definitions if td[0] == "global_options"]
    test_definitions = [td for td in test_definitions if td[0] != "global_options"]
    if len(global_options) == 0:
        global_options = ""
    elif len(global_options) == 1:
        global_options = global_options[0] + " "
    else:
        raise openloops.OpenLoopsError("Only one global_options specification is allowed.")
    override_options = " " + " ".join(args[2:])
else:
    global_options = []
    override_options = []



def run_process(test_def, dst_override=None, quiet=False):
    process_name = test_def[0]
    options = global_options + test_def[1] + override_options
    ko.parse(options)
    if dst_override:
        save = dst_override
    else:
        save = ko.save
    src, dst = data_files(None, save, process_name)
    if not dst and ko.save.lower() == "new":
        default_dst = default_data_file(process_name)
        if not os.path.exists(default_dst):
            dst = default_dst
    if not quiet:
        print "\nProcess:", process_name, "with options", options
        if ko.remaining_args:
            print "WARNING: unrecognised options:", ko.remaining_args
        print "Destination file:", dst
        print
    ptd = openloops.ProcessTestData(test_def[0],
                                    subprocess=ko.subprocess,
                                    phasespace=ko.energy,
                                    n=ko.n,
                                    parameters=openloops.Parameters(**ko.unknown_options),
                                    verbose=ko.verbose)
    if dst:
        ptd.dump(dst)
        print "Created reference data."



def check_process(test_def):
    process_name = test_def[0]
    options = global_options + test_def[1] + override_options
    print "\nProcess:", process_name, "with options", options, "(ignored in check mode)"
    ko.parse(options)
    if ko.remaining_args:
        print "WARNING: unrecognised options:", ko.remaining_args
    src, dst = data_files(ko.data, ko.save, process_name)
    src_msg = src
    if src and not os.path.exists(src):
        if ko.save.lower() == "new":
            dst = src
        src = None
        src_msg = "NOT FOUND"
    print "Source file:     ", src_msg
    print "Destination file:", dst
    print
    if src is None:
        run_process(test_def, dst_override=dst, quiet=True)
        if not dst:
            print "NOT VALIDATED"
    else:
        ptd = openloops.ProcessTestData.load(src)
        new_ptd = ptd.validate(verbose=ko.verbose)
        agrees = new_ptd.check_successful
        if agrees:
            print "Validation succeeded"
        else:
            print "VALIDATION FAILED"
        if dst is not None and (agrees or force):
            new_ptd.dump(dst)



for test_def in test_definitions:
    if check:
        check_process(test_def)
    else:
        run_process(test_def)
