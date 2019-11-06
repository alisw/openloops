
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
import atexit
import readline
import rlcompleter

import openloops_legacy as openloops

# command completion with tab
readline.parse_and_bind("tab: complete")

# command history
historyPath = os.path.expanduser("~/.openloops_history")

def save_history(historyPath=historyPath):
    import readline
    readline.write_history_file(historyPath)

if os.path.exists(historyPath):
    readline.read_history_file(historyPath)

atexit.register(save_history)

del os, atexit, readline, rlcompleter, save_history, historyPath

# change the prompt
sys.ps1 = "openloops> "
sys.ps2 = ".......... "

del sys

# attach openloops classes to the current namespace
Parameters = openloops.Parameters
PhaseSpacePoint = openloops.PhaseSpacePoint
Process = openloops.Process
MatrixElement = openloops.MatrixElement
ProcessTestData = openloops.ProcessTestData

# welcome
print """Welcome to the OpenLoops Python interface """ + ".".join(map(str,openloops.version)) + """.
This is a Python interpreter with the OpenLoops module loaded and library paths set."""
