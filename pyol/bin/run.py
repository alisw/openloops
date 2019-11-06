#!/usr/bin/env python
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

from __future__ import print_function

import sys
import argparse
import time
import openloops
import parallel

# =============== #
# argument parser #
# =============== #

def amptype_conv(a):
    a = a.lower()
    try:
        return int(a)
    except ValueError:
        return a

parser = argparse.ArgumentParser()
parser.add_argument('process', metavar='PROCESS',
                    help='The process to calculate')
parser.add_argument(
    '-a', '--amptype', type=amptype_conv, default=openloops.LOOP,
    choices=list(openloops.AMPTYPES.keys())+list(openloops.AMPTYPES.values()),
    help='amptype')
parser.add_argument(
    '-e', '--energy', type=float, default=openloops.default_energy,
    help='Energy to be used in phase space point generation.')
parser.add_argument(
    '-n', type=int, metavar='POINTS', default=None,
    help='Number of phase space points to calculate matrix elements for.')
parser.add_argument(
    '-t', '--time', dest='timing', action='store_true', default=False,
    help='''Measure runtime per phase space point
         (may be inaccurate for very simple processes).''')
parser.add_argument(
    '--tt', type=float, dest='mintime', metavar='MINTIME', default=10,
    help='Minimal time for runtime measurements in seconds.')
parser.add_argument(
    '--tn', type=int, dest='minn', metavar='MINN', default=1,
    help='Minimal number of points in runtime measurements.')
parser.add_argument(
    '-v', '--verbose', type=int, default=1, help='Verbosity level')
parser.add_argument(
    '-p', '--parallel', type=int, default=-1,
    help="""Number of subprocesses for parallel matrix element evaluation.
    Default -1 means no subprocesses, 0 means as many subprocesses as there
    are CPU cores. Note that 1 differs from -1 as it uses one subprocess
    instead of just the main process.""")
# This is just to create a help message, opt=val pairs are extracted manually
# to avoid ordering conflicts between options and positional arguments.
parser.add_argument(
    'options', metavar='OPT=VAL[:VAL2[...]]', nargs='*',
    help="""Options to pass to directly to OpenLoops. Swap options:
    separate multiple values for the same option with colons to evaluate the
    matrix element for the same phase space point with different options.
    The first value of each swap option is initialised before registering
    processes, i.e. the first values should be such that the loaded process
    works with all following values. Can be used with multiple options if they
    all have the same number of swap options.""")


def set_option(key, val):
    if key.startswith('alpha') and '/' in val:
        try:
            valnum, valden = val.split('/')
            val = float(valnum)/float(valden)
        except ValueError:
            print('[PYOL] ERROR: invalid option \'{}\''.format(opt))
    if args.verbose >= 3:
        print('call set_parameter({},{})'.format(key,val))
    openloops.set_parameter(key, val)


def eval_me(proc, psp):
    mes = []
    if swap_options:
        # multiple evaluation with different parameters
        for sopts in swap_options:
            for key, val in sopts:
                set_option(key, val)
            mes.append(proc.evaluate(psp))
    else:
        # single evaluation
        mes.append(proc.evaluate(psp))
    return mes


def print_me(mes, first=False):
    """Print matrix elements.
    mes -- a list of matrix elements,
           usually for the same phase space point.
    first -- if True, print a line to lable the results.
    """
    if first and args.verbose >= 1:
        prepend = ''
        if swap_options:
            prepend = '   '
        if mes[0].amptype == openloops.TREE:
            print(prepend + '{:>23}'.format('tree'))
        elif mes[0].amptype == openloops.LOOP:
            print((prepend + 5*'{:>23}').format(
                  'tree', 'finite', 'ir1', 'ir2', 'acc'))
        elif mes[0].amptype == openloops.LOOP2:
            print((prepend + 2*'{:>23}').format('loop2', 'acc'))
    if args.verbose >= 2:
        print(mes[0].psp)
    if args.verbose >= 1:
        if swap_options:
            for ns, me in enumerate(mes):
                print('[{}]{}'.format(ns+1, me.valuestr()))
        else:
            print(mes[0].valuestr())

# =============== #
# option handling #
# =============== #

stroptions = openloops.config['pyrunoptions']
stroptions.extend([arg for arg in sys.argv[1:]
                 if ('=' in arg and not arg.startswith('-'))])
args = parser.parse_args(
    [arg for arg in sys.argv[1:] if (arg.startswith('-') or '=' not in arg)])
if not args.timing and args.n is None:
    args.n = 1

options = [] # [(opt, val), ...]
swap_options = [] # [(opt, [val1, val2, ...]), ...]

for opt in stroptions:
    try:
        key, val = opt.split('=',1)
    except ValueError:
        print('[PYOL] ERROR: invalid option \'{}\''.format(opt))
        sys.exit(1)
    vals = val.split(':')
    options.append((key,vals[0]))
    if len(vals) > 1:
        if not swap_options:
            swap_options = [list() for val in vals]
        if len(vals) != len(swap_options):
            print('ERROR: swap options must all have equal length.')
            sys.exit(1)
        for val, sopt in zip(vals, swap_options):
            sopt.append((key,val))

# ======================== #
# parameter initialisation #
# ======================== #

# Also in parallel mode parameters must be initialised in the main process,
# so that the processes can be registered to generate phase space points
# in the main programm (which are then distributed to the workers).

for key, val in options:
    set_option(key, val)

# ================== #
# register processes #
# ================== #

# Note that when the processes are successfully registered in the main process
# we assume that this also works in the subprocesses.

is_library = False
if not '>' in args.process:
    is_library = True
    try:
        procinfo = openloops.ProcessInfo(args.process)
    except openloops.ProcessInfoError:
        is_library = False

if not is_library:
    # register a partonic channel
    try:
        processes = [openloops.Process(args.process, args.amptype)]
    except openloops.RegisterProcessError:
        print(('[PYOL] ERROR registering process \'{}\' ' +
               'with amptype {}.').format(args.process, args.amptype))
        if not '>' in args.process:
            print('(reading info file failed, too)')
        sys.exit(1)
else:
    # register channels from a process library
    try:
        processes = procinfo.register()
    except openloops.RegisterProcessError as err:
        print('[PYOL] ERROR while registering process from library ' +
              '{}:\n       {}').format(args.process, err.args[0])
        sys.exit(1)

# === #
# run #
# === #

if args.parallel < 0:
    # evaluation in main process only
    for proc in processes:
        print()
        print('"' + proc.process + '"')
        if not args.timing:
            for n in range(args.n):
                psp = proc.psp(args.energy)
                mes = eval_me(proc, psp)
                print_me(mes, first=not(n))
        else:
            # timing
            psp = proc.psp(args.energy)
            mes = eval_me(proc, psp)
            print_me(mes, first=True)
            psp = proc.psp(args.energy)
            mes = eval_me(proc, psp)
            print_me(mes)
            starttime = time.clock()
            npoints = 0
            if args.n is not None:
                npoints = args.n
                for n in range(args.n):
                    psp = proc.psp(args.energy)
                    mes = eval_me(proc, psp)
                    print_me(mes)
            else:
                while (time.clock() < starttime + args.mintime or
                    npoints < args.minn):
                    npoints = npoints + 1
                    psp = proc.psp(args.energy)
                    mes = eval_me(proc, psp)
                    print_me(mes)
            endtime = time.clock()
            print(('time per phase space point: {:3f} ms (avg. of {} ' +
                'points)').format(1000*(endtime - starttime)/npoints, npoints))

else:
    # parallel evaluation in subprocesses
    procamp_list = [(proc.process, proc.amptype) for proc in processes]
    pool = parallel.Pool(procamp_list, args.parallel, options, swap_options)
    for proc in processes:
        print()
        print('"' + proc.process + '"')
        # submit phase space points
        for n in range(args.n):
            pool.put((proc.process, proc.amptype), proc.psp(args.energy))
        for n in range(args.n):
            mes = pool.get()
            print_me(mes, first=not(n))
