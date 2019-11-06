
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


ko = keywordoptions.KeywordOptions()

ko.add("subprocess", converter=openloops.parse_subprocess_argument)
ko.add("verbose", converter=int, default=1)

ko.parse(sys.argv[2:])

data_file = sys.argv[1]

try:
    ptd = openloops.ProcessTestData.load(data_file)
except IOError:
    ptd = openloops.ProcessTestData.load(os.path.join(openloops.test_data_folder, data_file + ".ptd"))

print
ptd.show(verbose=ko.verbose)
