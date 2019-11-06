#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#!******************************************************************************!
#! Copyright (C) 2014-2019 OpenLoops Collaboration. For authors see authors.txt !
#!                                                                              !
#! This file is part of OpenLoops.                                              !
#!                                                                              !
#! OpenLoops is free software: you can redistribute it and/or modify            !
#! it under the terms of the GNU General Public License as published by         !
#! the Free Software Foundation, either version 3 of the License, or            !
#! (at your option) any later version.                                          !
#!                                                                              !
#! OpenLoops is distributed in the hope that it will be useful,                 !
#! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
#! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
#! GNU General Public License for more details.                                 !
#!                                                                              !
#! You should have received a copy of the GNU General Public License            !
#! along with OpenLoops.  If not, see <http://www.gnu.org/licenses/>.           !
#!******************************************************************************!

# This module can be run as a command line tool
# or imported as a module for interactive use.

# TODO
# * Show more trigger information in 2x histograms.

# http://stackoverflow.com/questions/7534453/matplotlib-does-not-show-my-drawings-although-i-call-pyplot-show

from __future__ import division
import os
import argparse
import math
import collections
import tarfile
try:
    from matplotlib import pyplot
    from matplotlib.ticker import MultipleLocator
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("""*** matplotlib is not available. ***
Unless matplotlib is installed, creating plots is not possible.
Only data accumulation can be performed.""")

try:
    strtype = basestring
except NameError:
    strtype = str

# lower limit of the x axis
xlimit_lower = -16
# minimal upper limit of the x axis; will be extended if more data is available
xlimit_upper = 3


class StabilityData(object):

    def __init__(self, total=0, total_qp=0, points=None, points_qp=None,
                 qp=0, killed=0, killed_qp=0, channel=None):
        # Use separate values for dp and qp total/killed points, because
        # the numbers are both in dp and qp files and one might want to
        # import only dp or qp files.
        # total_qp might become useful in 2x-mode histograms.
        self.total = total
        self.total_qp = total_qp
        self.points = points
        self.points_qp = points_qp
        self.qp = qp
        self.killed = killed
        self.killed_qp = killed_qp
        self.channel = channel

    def add(self, total=0, total_qp=0, points=None, points_qp=None,
            qp=0, killed=0, killed_qp=0, channel=None):
        """Add histogram data."""
        if points:
            self.total += total
            if self.points is None:
                self.points = points
            else:
                self.points = list(map(sum, zip(self.points, points)))
            self.killed += killed
        if points_qp:
            self.total_qp += total_qp
            if self.points_qp is None:
                self.points_qp = points_qp
            else:
                self.points_qp = list(map(sum, zip(self.points_qp, points_qp)))
            self.qp += qp
            self.killed_qp += killed_qp
        if channel:
            if not self.channel:
                self.channel = channel
            elif self.channel != channel:
                raise Exception('StabilityData.add(): cannot add data ' +
                                'from different channels')

    def add_dataline(self, dataline, qp=False, channel=None):
        """Add histogram data from a string.
        Format: 'total | h1 h2 h3 ... | qp killed'"""
        dat = dataline.split('|')
        try:
            total = int(dat[0])
            points = list(map(int, dat[1].split()))
            triggerdat = list(map(int, dat[2].split()))
            if len(triggerdat) == 8:
                # workaround for missing np(2:6) slicing
                # in old finish_histograms()
                triggerdat[1:6]
            n_qp, killed = triggerdat[-2:]
        except:
            print('unexpected histogram data format:')
            print('\'' + dataline + '\'')
            raise
        if qp:
            self.add(None, total, None, points, n_qp, None, killed, channel)
        else:
            self.add(total, None, points, None, n_qp, killed, None, channel)

    def export(self, outdir='', channel=None):
        """Export dp and qp stability data to files
        which are readable by import_files(). Options:
        outdir -- output directory, default: working directory
        channel -- channel name, needed only if it was not extracted
                   from the names of imported files."""
        if channel is None:
            if self.channel is None:
                raise Exception('StabilityData.export(): channel name needed')
            else:
                channel = self.channel
        filename = os.path.join(outdir, 'histogram_' + channel + '_acc.log')
        filename_qp = os.path.join(outdir,
                                   'histogram_' + channel + '_qp_acc.log')
        with open(filename, 'w') as fh:
            fh.write(str(self.total))
            fh.write(' | ')
            fh.write(' '.join(map(str, self.points)))
            fh.write(' | ')
            fh.write(' '.join(map(str, [self.qp, self.killed])))
            fh.write('\n')
        with open(filename_qp, 'w') as fh:
            fh.write(str(self.total_qp))
            fh.write(' | ')
            fh.write(' '.join(map(str, self.points_qp)))
            fh.write(' | ')
            fh.write(' '.join(map(str, [self.qp, self.killed_qp])))
            fh.write('\n')

    def __add__(self, other):
        if self.channel != other.channel:
            raise Exception('StabilityData: only equal channels can be summed')
        if self.channel:
            channel = self.channel
        else:
            channel = other.channel
        return StabilityData(total = self.total + other.total,
                             total_qp = self.total_qp + other.total_qp,
                             points = self.points + other.points,
                             points_qp = self.points_qp + other.points_qp,
                             qp = self.qp + other.qp,
                             killed = self.killed + other.killed,
                             killed_qp = self.killed_qp + other.killed_qp,
                             channel = channel)

    def __repr__(self):
        return ('StabilityData(total={}, total_qp={}, points={}, ' +
                'points_qp={}, qp={}, killed={}, killed_qp={}, channel=\'{}\')'
                ).format(self.total, self.total_qp,
                         self.points, self.points_qp, self.qp,
                         self.killed, self.killed_qp, self.channel)


def select_file(filename, libraries, channels):
    """If filename is a valid histogram data file name,
    return the corresponding library name, channel name and
    True (False) if it does (not) contain quad precision data.
    Otherwise return (None, None, None) for invalid file names
    and (False, False, False) if the file is not selected."""
    f = os.path.basename(filename)
    if f.startswith('histogram_') and f.endswith('.log'):
        splitf = os.path.basename(f).split('_')
        if splitf[-2] == 'qp':
            lib = '_'.join(splitf[1:-4])
            ch = '_'.join(splitf[-4:-2])
            qp = True
        else:
            lib = '_'.join(splitf[1:-3])
            ch = '_'.join(splitf[-3:-1])
            qp = False
        if ((libraries and lib not in libraries) or
              (channels and ch not in channels)):
            lib = False
            ch = False
            qp = False
        else:
            ch = lib + '_' + ch
        return (lib, ch, qp)
    else:
        return (None, None, None)


def import_files(filesdirs, libraries=None, channels=None):
    """Import stability data from files, directories and tar archives.
    Return a list of StabilityData objects for each channel.
    Only the first line in each file is considered.
    Options:
    libraries -- list of libraries to select; None or empty for all.
    channels -- list of channels to select; None or empty for all."""
    if isinstance(filesdirs, strtype):
        filesdirs = [filesdirs]
    if channels:
        channels = [ch + '_1' if len(ch.split('_')) == 1 else ch
                    for ch in channels]
    data_dict = collections.defaultdict(StabilityData)

    for filedir in filesdirs:
        # add all files from the directory
        if os.path.isdir(filedir):
            no_imported_file = True
            for f in os.listdir(filedir):
                lib, ch, qp = select_file(f, libraries, channels)
                if not lib:
                    continue
                no_imported_file = False
                with open(os.path.join(filedir, f)) as fh:
                    data_dict[ch].add_dataline(
                        fh.readline(), qp=qp, channel=ch)
            if no_imported_file:
                print("""Directory '{}'
does not contain any histogram files which are suited
for automatic channel name extraction.""".format(f))
        # add all files from the tar archive
        elif os.path.isfile(filedir) and (filedir.endswith('.tar') or
             filedir.endswith('.tar.gz') or filedir.endswith('.tar.bz2') or
             filedir.endswith('.tgz')):
            with tarfile.open(filedir) as tar:
                no_imported_file = True
                for f in tar:
                    if f.isfile():
                        lib, ch, qp = select_file(f.name, libraries, channels)
                        if not lib:
                            continue
                        no_imported_file = False
                        tarfh = tar.extractfile(f)
                        data_dict[ch].add_dataline(
                            tarfh.readline(), qp=qp, channel=ch)
                        tarfh.close()
                if no_imported_file:
                    print("""Archive '{}'
does not contain any histogram files which are suited
for automatic channel name extraction.""".format(f))
        # add individual files
        elif os.path.isfile(filedir):
            lib, ch, qp = select_file(filedir, libraries, channels)
            if lib is None:
                raise Exception('File name unsuited for automatic channel ' +
                                'name extraction: \'{}\''.format(filedir))
            if not lib:
                continue
            with open(filedir) as fh:
                data_dict[ch].add_dataline(
                    fh.readline(), qp=qp, channel=ch)

        else:
            raise Exception(
                '\'{}\' is neither a file nor a directory'.format(filedir))

    return [data_dict[key] for key in sorted(data_dict.keys())]


def fraction_format(x):
    if round(100*x) >= 100:
        return '{}%'.format(round(100*x))
    elif x < 0.0001:
        return '{:.2g}ppm'.format(x*1000000)
    elif x < 0.001:
        return '{:.3g}ppm'.format(x*1000000)
    else:
        return '{:.2g}%'.format(x*100)


def stability_plot(data):
    """Add a step plot with data from a StabilityData object."""
    pyplot.figure()
    # usually 19 data points, starting with the bin between -16 and -15.
    # The points in the bin between 1 and 2 deviate by a factor 10 to 100.
    if data.channel:
        pyplot.suptitle(data.channel, fontweight='bold')
    killed = max(data.killed, data.killed_qp)
    if data.total_qp:
        killed = data.killed_qp
    else:
        killed = data.killed
    if data.total:
        total = data.total
    else:
        total = data.total_qp
    pyplot.title('points: {}, qp: {} ({}), killed: {} ({})'.format(
        total,
        data.qp, fraction_format(data.qp/total),
        killed, fraction_format(killed/total)),
        fontsize=10)
    pyplot.xlabel('estimated precision')
    pyplot.ylabel('fraction of points, cumulative')
    # double precision histogram
    if data.points:
        dp_points = list(data.points)
        while dp_points and dp_points[-1] == 0:
            dp_points.pop()
        if not dp_points:
            print('Note: histogram for channel ' + data.channel + 'is empty. Skipped.')
            pyplot.clf()
            pyplot.close()
            return False
        # If a value is zero in a logarithmic step plot, matplotlib
        # will omit the vertical line --> use a small number instead.
        n_dp_points = len(dp_points)
        # minimal non-zero data point
        min_dp_point = min(y for y in dp_points if y > 0)/total
        # append zero to make the last data point visible
        dp_points.append(0)
        # If a value is zero in a logarithmic step plot, matplotlib
        # will omit the vertical line --> use a small number instead.
        dp_points = [n if n > 0 else 1e-99 for n in dp_points]
        pyplot.step(range(-16,len(dp_points)-16),
                    [y/total for y in dp_points], where='post')
    else:
        n_dp_points = 0
        min_dp_point = None
    # the same for the quad precision histogram
    if data.points_qp and data.points_qp[0]:
        qp_points = list(data.points_qp)
        while qp_points[-1] == 0:
            qp_points.pop()
        n_qp_points = len(qp_points)
        min_qp_point = min(y for y in qp_points if y > 0)/total
        qp_points.append(0)
        qp_points = [n if n > 0 else 1e-99 for n in qp_points]
        pyplot.step(range(-16,len(qp_points)-16),
                    [y/total for y in qp_points], where='post')
    else:
        n_qp_points = 0
        min_qp_point = None

    if min_dp_point is None:
        min_point = min_qp_point
    elif min_qp_point is None:
        min_point = min_dp_point
    else:
        min_point = min(min_dp_point, min_qp_point)
    xlimit_upper_local = max(
        xlimit_upper,
        n_dp_points + xlimit_lower, n_qp_points + xlimit_lower)
    pyplot.xlim(xlimit_lower, xlimit_upper_local)
    pyplot.xticks(range(xlimit_lower, xlimit_upper_local, 2))
    # lower limit of y axis values for log scale:
    # 10^n with the largest integer n such that 10^n is smaller
    # than the smallest value
    miny = 10**math.floor(math.log(min_point, 10))
    pyplot.ylim(miny, 2)
    # minor ticks at integer x values
    pyplot.gca().xaxis.set_minor_locator(MultipleLocator(1))
    pyplot.semilogy()
    return True


if __name__ == '__main__':
    # option handling
    parser = argparse.ArgumentParser(
        description="""Create stability histogram plots from OpenLoops
stability log files. Files can be specified directly, as (first-level) content
of directories or as tar archives. File names must follow the OpenLoops naming
convention. If directories or archives are included, files which do not match
the naming convention are skipped silently. Plots are shown unless --output is
specified. Showing plots only works if the matplotlib back end is properly
configured and the back end is available from the terminal which is used.""")
    parser.add_argument('-o', '--output', dest='outfile', metavar='outfile',
                        default='', help='PDF output file name.')
    parser.add_argument('-a', '--accumulate', action='store_true', default=False,
                        help="""Save accumulated histogram data to files instead
of creating plots. Per channel: histogram_<channel>_[qp_]acc.log. The output
directory can be specified with --output (directory must exist).""")
    parser.add_argument('-c', '--channel', action='append', dest='channels',
                        metavar='channel',
                        help="""Select the specified partonic channel from the
input files. Can be specified multiple times.""")
    parser.add_argument('-l', '--library', action='append', dest='libraries',
                        metavar='library',
                        help="""Select files from the specified library only.
Can be specified multiple times.""")
    parser.add_argument('files', metavar='file_or_directory', nargs='+',
                        help='File or directory to import.')
    args = parser.parse_args()

    outfile = args.outfile
    if outfile:
        if not outfile.lower().endswith('.pdf'):
            outfile += '.pdf'

    # read histogram data and create plots
    stability_data = import_files(args.files, libraries=args.libraries,
                                  channels=args.channels)
    if args.accumulate:
        for data in stability_data:
            data.export(outdir=args.outfile)
    else:
        plots_created = 0
        for data in stability_data:
            if stability_plot(data):
                plots_created = plots_created + 1
        # output
        if outfile:
            pdf = PdfPages(outfile)
            for n in range(plots_created):
                pyplot.figure(n+1)
                pdf.savefig()
            pdf.close()
        else:
            pyplot.show()
