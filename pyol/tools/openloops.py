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


# TODO
# * ProcessInfo(): select loop/tree processes only,
#                  use Type= to group channels.
# * If config is loaded, replace loopspec_flags.
# * Retrieve install_path (needs getparameter_string)?
# * Take proclib_dir from OL config?
# * proclib_dir is currently relative to the working directory.

from __future__ import print_function

import os
import sys
import atexit
import ctypes
import collections
from ctypes import c_int, c_double, c_char_p, byref, POINTER
import OLBaseConfig

try:
    strtype = basestring
except NameError:
    strtype = str

config = OLBaseConfig.get_config()

if config['print_python_version']:
    print('openloops.py uses Python', sys.version)

proclib_dir = 'proclib'

c_int_ptr = POINTER(c_int)
c_double_ptr = POINTER(c_double)

TREE = 1
CC = 2
SC = 3
SCPV = 4
LOOP = 11
LOOP2 = 12
AMPTYPES = {'tree': TREE, 'cc': CC, 'sc': SC, 'scpv': SCPV,
            'loop': LOOP, 'loop2': LOOP2}
loopspec_flags = 'tlsp'

default_energy = 1000
default_amptype = LOOP
numberformat = '{:>23.15e}'

class OpenLoopsError(Exception):
    pass

class RegisterProcessError(OpenLoopsError):
    pass

class ProcessInfoError(OpenLoopsError):
    pass

#from ctypes.util import find_library
#libopenloops = find_library('openloops')
#if not libopenloops:
    #print('ERROR: openloops library not found')
    #sys.exit(1)
#ol = ctypes.CDLL(libopenloops)

if sys.platform.startswith("darwin"):
    ol = ctypes.CDLL("libopenloops.dylib")
else:
    ol = ctypes.CDLL("libopenloops.so")

class LibraryContent(object):
    def __init__(self, content, func=None, args=None):
        # content bit flag order: 1=tree 2=loop 4=loops 8=pt
        self.tree = bool(content & 1)
        self.loop = bool((content >> 1) & 1)
        self.loop2 = bool((content >> 2) & 1)
        self.pt = bool((content >> 3) & 1)

# (const char* key, int val)
setparameter_int_c = ol.ol_setparameter_int
setparameter_int_c.argtypes = [c_char_p, c_int]
setparameter_int_c.restype = None
def setparameter_int(key, val):
    setparameter_int_c(c_char_p(key.encode()), val)

setparameter_bool_c = ol.ol_setparameter_bool
setparameter_bool_c.argtypes = [c_char_p, c_int]
setparameter_bool_c.restype = None
def setparameter_bool(key, val):
  if val:
    setparameter_bool_c(c_char_p(key.encode()), 1)
  else:
    setparameter_bool_c(c_char_p(key.encode()), 0)

# (const char* key, double val)
setparameter_double_c = ol.ol_setparameter_double
setparameter_double_c.argtypes = [c_char_p, c_double]
setparameter_double_c.restype = None
def setparameter_double(key, val):
    setparameter_double_c(c_char_p(key.encode()), val)

# (const char* key, char* val)
setparameter_string_c = ol.ol_setparameter_string
setparameter_string_c.argtypes = [c_char_p, c_char_p]
setparameter_string_c.restype = None
def setparameter_string(key, val):
    setparameter_string_c(c_char_p(key.encode()), c_char_p(val.encode()))

# (const char* key, int[1] val)
getparameter_int_c = ol.ol_getparameter_int
getparameter_int_c.argtypes = [c_char_p, c_int_ptr]
getparameter_int_c.restype = None

# (const char* key, double[1] val)
getparameter_double_c = ol.ol_getparameter_double
getparameter_double_c.argtypes = [c_char_p, c_double_ptr]
getparameter_double_c.restype = None

# (const char* proc, int amptype) -> int
register_process_c = ol.ol_register_process
register_process_c.argtypes = [c_char_p, c_int]
register_process_c.restype = int
def register_process(proc, amptype):
    return register_process_c(c_char_p(proc.encode()), amptype)

# (int id) -> int
n_external = ol.ol_n_external
n_external.argtypes = [c_int]
n_external.restype = int

# (int id) -> int
amplitudetype = ol.ol_amplitudetype
n_external.argtypes = [c_int]
n_external.restype = int

# (int id, double sqrt_s, double[5*n] pp)
phase_space_point_c = ol.ol_phase_space_point
phase_space_point_c.argtypes = [c_int, c_double, c_double_ptr]
phase_space_point_c.restype = None

parameters_flush = ol.ol_parameters_flush
parameters_flush.argtypes = []
parameters_flush.restype = None

parameters_write = ol.ol_printparameter
parameters_write.argtypes = []
parameters_write.restype = None

start = ol.ol_start
start.argtypes = []
start.restype = None
start.started = False

finish = ol.ol_finish
finish.argtypes = []
finish.restype = None

library_content = ol.ol_library_content
library_content.argtypes = [c_int]
library_content.restype = int
library_content.errcheck = LibraryContent

# (int id, double[5*n] pp, double[1] tree)
evaluate_tree_c = ol.ol_evaluate_tree
evaluate_tree_c.argtypes = [c_int, c_double_ptr, c_double_ptr]
evaluate_tree_c.restype = None

# (int id, double[5*n] pp, double[1] tree, double[r] cc, double[1] ewcc)
evaluate_cc_c = ol.ol_evaluate_cc
evaluate_cc_c.argtypes = [c_int, c_double_ptr, c_double_ptr,
                          c_double_ptr, c_double_ptr]
evaluate_cc_c.restype = None

# (int id, double[5*n] pp, int emitter, double[4] polvect, double[n] sc)
evaluate_sc_c = ol.ol_evaluate_sc
evaluate_sc_c.argtypes = [
    c_int, c_double_ptr, c_int, c_double_ptr, c_double_ptr]
evaluate_sc_c.restype = None

# (int id, double[5*n] pp, double[1] tree, double[3] loop, double[3] ir1,
#  double[5] loop2, double[5] ir2, double[1] acc)
evaluate_full_c = ol.ol_evaluate_full
evaluate_full_c.argtypes = [
    c_int, c_double_ptr, c_double_ptr, c_double_ptr,
    c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr]
evaluate_full_c.restype = None

# (int id, double[5*n] pp, double[1] tree, double[3] loop, double[1] acc)
evaluate_loop_c = ol.ol_evaluate_loop
evaluate_loop_c.argtypes = [c_int, c_double_ptr, c_double_ptr,
                            c_double_ptr, c_double_ptr]
evaluate_loop_c.restype = None

# (int id, double[5*n] pp, double[1] loop2, double[1] acc)
evaluate_loop2_c = ol.ol_evaluate_loop2
evaluate_loop2_c.argtypes = [
    c_int, c_double_ptr, c_double_ptr, c_double_ptr]
evaluate_loop2_c.restype = None

# (int id, double[5*n] pp, double[1] tree, double[1] ct)
evaluate_ct_c = ol.ol_evaluate_ct
evaluate_ct_c.argtypes = [
    c_int, c_double_ptr, c_double_ptr, c_double_ptr]
evaluate_ct_c.restype = None

# (int id, double[5*n] pp, double[1] tree, double[1] r2)
evaluate_r2_c = ol.ol_evaluate_r2
evaluate_r2_c.argtypes = [
    c_int, c_double_ptr, c_double_ptr, c_double_ptr]
evaluate_r2_c.restype = None

# (int id, double[5*n] pp, double[1] tree, double[1] pt, double[1] loop)
evaluate_pt = ol.ol_evaluate_pt
evaluate_pt.argtypes = [c_int, c_double_ptr, c_double_ptr,
                        c_double_ptr, c_double_ptr]
evaluate_pt.restype = None

# ol_tree_colbasis_dim(int id, int[1] ncolb, int[1] colelemsz, int[1] nhel)
# ol_tree_colbasis(int id, int* basis, int* needed)
#   basis(tree_colbasis_elemsize(id),get_tree_colbasis_dim(id))
#   needed(get_tree_colbasis_dim(id),get_tree_colbasis_dim(id))
# ol_evaluate_tree_colvect(int id, double[5*n] pp, double* amp, int[1] nhel)
#   amp(2*get_tree_colbasis_dim(id),get_nhel(id))
# ol_evaluate_ccmatrix(int id, double[5*n] pp, double[1] tree, double[n][n] ccij, double[1] ewcc)
# ol_evaluate_scpowheg(int id, double[5*n] pp, int emitter, double[1] res, double[n][n] resmunu)

atexit.register(finish)


def set_parameter(key, val):
    if type(val) == int:
        setparameter_int(key, val)
    elif type(val) == bool:
        setparameter_bool(key, val)
    elif isinstance(val, float):
        setparameter_double(key, val)
    elif isinstance(val, strtype):
        if key.startswith('alpha') and '/' in val:
            try:
                valnum, valden = val.split('/')
                val = float(valnum) / float(valden)
            except ValueError:
                raise OpenLoopsError(
                    'Invalid option \'{}={}\''.format(key, val))
            setparameter_double(key, val)
        else:
            setparameter_string(key, val)
    else:
        raise TypeError('set_parameter() value argument must be int, ' +
                        'float or str, not \'{}\''.format(type(val)))

def get_parameter_int(key):
    val_c = c_int()
    getparameter_int_c(c_char_p(key.encode()), byref(val_c))
    return val_c.value

def get_parameter_double(key):
    val_c = c_double()
    getparameter_double_c(c_char_p(key.encode()), byref(val_c))
    return val_c.value


class PhaseSpacePoint(object):
    """A phase space point, stored as a 5*n c_double array.
    Can be passed to C function handles directly.
    The 5th momentum component can store the mass, but is neither used by
    OpenLoops nor calculated from the momentum by this class."""
    def __init__(self, pp, n=None):
        """PhaseSpacePoint constructor. Valid input formats for n particles:
        * 5*n c_double array;
        * 5*n list;
        * 4*n list if n is given;
        * 2-dimenional list of 4 (or 5) component momenta.
        Note: does not generate random points."""
        if isinstance(pp, (list, tuple)):
            if isinstance(pp[0], (list, tuple)) or n:
                if not isinstance(pp[0], (list, tuple)):
                    # 1-dimensional -> 2-dimensional using n
                    if len(pp) % n:
                        raise OpenLoopsError('Invalid phase space point.')
                    k = len(pp)/n
                    pp = [(pp[m:m+k]) for m in range(0,len(pp),k)]
                ppls = []
                if n and len(pp) != n:
                    raise OpenLoopsError(
                        'Invalid phase space point: wrong number of momenta')
                for p in pp:
                    ppls.extend(p)
                    if len(p) == 4:
                        ppls.append(-1)
                    elif len(p) != 5:
                        raise OpenLoopsError(
                            'Invalid phase space point: momentum must have 4' +
                            ' or 5 components.')
            else:
                # 1-dimensional and n not known: assume 5*n notation
                ppls = pp
        else:
            # assume that a valid c_double array was passed;
            # make a copy, otherwise the data can be overwritten
            ppls = tuple(pp)
        self.pp = (len(ppls)*c_double)(*ppls)
        self._as_parameter_ = ctypes.cast(self.pp,c_double_ptr)
    def __getstate__(self):
        # pickle as tuple
        return tuple(self.pp)
    def __setstate__(self, ppls):
        self.pp = (len(ppls)*c_double)(*ppls)
        self._as_parameter_ = ctypes.cast(self.pp,c_double_ptr)
    def __str__(self):
        return str(tuple(self.pp))


# The variable name must be the same as the name of the namedtuple,
# otherwise it is not pickleable.
LoopME = collections.namedtuple('LoopME', ['finite', 'ir1', 'ir2'])
IOperator = collections.namedtuple('IOperator', ['finite', 'ir1', 'ir2'])
Loop2ME = collections.namedtuple('Loop2ME',
                                 ['finite', 'ir1', 'ir2', 'ir3', 'ir4'])


class MatrixElement(object):
    def __init__(self, amptype, psp, **kwargs):
        self.amptype = amptype
        self.psp = psp
        for key, val in kwargs.items():
            setattr(self, key, val)
    def __str__(self):
        if self.amptype == TREE:
            return 'tree={}'.format(self.tree)
        elif self.amptype == LOOP:
            return 'tree={} {} acc={}'.format(self.tree, self.loop, self.acc)
        elif self.amptype == LOOP2:
            return 'loop2={} acc={}'.format(self.loop2.finite, self.acc)
        elif self.amptype == 'ct':
            return 'tree={} ct={}'.format(self.tree, self.ct)
        elif self.amptype == 'r2':
            return 'tree={} r2={}'.format(self.tree, self.r2)
        else:
            return ('MatrixElement.__str__ not implemented for amptype {}'
                   ).format(self.amptype)
    def valuestr(self):
        if self.amptype == TREE:
            return numberformat.format(self.tree)
        elif self.amptype == LOOP:
            return (5*numberformat).format(self.tree, self.loop.finite,
                                           self.loop.ir1, self.loop.ir2,
                                           self.acc)
        elif self.amptype == LOOP2:
            return (2*numberformat).format(self.loop2.finite, self.acc)
        elif self.amptype == 'ct':
            return (2*numberformat).format(self.tree, self.ct)
        elif self.amptype == 'r2':
            return (2*numberformat).format(self.tree, self.r2)
        else:
            return ('MatrixElement.valuestr not implemented for amptype {}'
                   ).format(self.amptype)


class Process(object):

    _tree_buf = c_double()
    _ewcc_buf = c_double()
    _ct_buf = c_double()
    _pt_buf = c_double()
    _ptloop_buf = c_double()
    _loop_buf = (3*c_double)()
    _ir1_buf = (3*c_double)()
    _loop2_buf = (5*c_double)()
    _ir2_buf = (5*c_double)()
    _acc_buf = c_double()

    def __init__(self, process, amptype=default_amptype, qcd=-1, ew=-1):
        """Register a process with given amptype (default=LOOP).
        Use 'qcd' rsp. 'ew' to request a process in a certain order
        in alpha_s rsp. alpha.

        public members:
        * process -- the process name as registered
        * id (int) -- process id
        * amptype (int) -- amptype
        * contains -- LibraryContent() instance
        * n -- number of external particles
        Process instances can be passed to C functions as id.
        """
        typearg = amptype
        needs_tree = False
        needs_pt = False
        if isinstance(amptype, strtype):
            try:
                amptype = AMPTYPES[amptype.lower()]
            except KeyError:
                pass
        if isinstance(amptype, strtype):
            if not set(typearg) - set(loopspec_flags):
                if 's' in typearg and 'l' in typearg:
                    amptype = LOOP2
                elif 'l' in amptype:
                    amptype = LOOP
                elif 't' in amptype:
                    amptype = TREE
            if 't' in typearg:
                needs_tree = True
            if 'p' in typearg:
                needs_pt = True
        if not amptype in AMPTYPES.values():
            raise RegisterProcessError(
                'Process(): illegal amptype \'{}\''.format(typearg))
        if qcd >= 0:
            set_parameter('order_qcd', qcd)
        if ew >= 0:
            set_parameter('order_ew', ew)
        self.id = register_process(process, amptype)
        if self.id <= 0:
            raise RegisterProcessError(
                'Failed to load process \'{}\' with amptype {}'.format(
                process, typearg))
        self.contains = library_content(self.id)
        if needs_tree and not self.contains.tree:
            raise RegisterProcessError((
                'Dedicated tree matrix elements are not available ' +
                'in process {}.').format(process))
        if needs_pt and not self.contains.pt:
            raise RegisterProcessError((
                'Pseudo-tree matrix elements are not available ' +
                'in process {}.').format(process))
        self.process = process
        self.amptype = amptype
        self.n = n_external(self.id)
        self._as_parameter_ = self.id
        self._pp_buf = ((5*self.n) * c_double)()
        if amptype == CC:
            self._cc_buf = (((self.n*(self.n-1))/2) * c_double)()
        if amptype == SC:
            self._sc_buf = ((2*self.n**2) * c_double)()
        if amptype == SCPV:
            self._sc_polvect_buf = ((self.n) * c_double)()

    def psp(self, sqrt_s=default_energy):
        """Generate a random phase space point for the process."""
        phase_space_point_c(self.id, sqrt_s, self._pp_buf)
        return PhaseSpacePoint(self._pp_buf, self.n)

    def evaluate(self, pp_or_sqrt_s=default_energy, amptype=None):
        """Calculate matrix elements either for a given phase space point
        or for a random phase space point of given or default energy.

        'amptype' overrides the amptype of the process -- if you use this,
                  it's up to you to figure out if what you do makes sense.
                  'ct' calculates the counterterms only."""
        ct_only = False
        r2_only = False
        if amptype == 'ct':
          ct_only = True
          amptype = None
        elif amptype == 'r2':
          r2_only = True
          amptype = None
        if not amptype:
            amptype = self.amptype
        if not start.started:
            start()
            start.started = True
        if isinstance(pp_or_sqrt_s, (int, float)):
            psp = self.psp(pp_or_sqrt_s)
        elif isinstance(pp_or_sqrt_s, PhaseSpacePoint):
            psp = pp_or_sqrt_s
        else:
            psp = PhaseSpacePoint(pp_or_sqrt_s, self.n)

        if amptype in (LOOP, LOOP2):
            if ct_only:
                evaluate_ct_c(self.id, psp, Process._tree_buf, Process._ct_buf)
                me = MatrixElement('ct', psp, tree=Process._tree_buf.value,
                                   ct=Process._ct_buf.value)
            elif r2_only:
                evaluate_r2_c(self.id, psp, Process._tree_buf, Process._ct_buf)
                me = MatrixElement('r2', psp, tree=Process._tree_buf.value,
                                   r2=Process._ct_buf.value)
            else:
                evaluate_full_c(self.id, psp, Process._tree_buf,
                                Process._loop_buf, Process._ir1_buf,
                                Process._loop2_buf, Process._ir2_buf,
                                Process._acc_buf)
                me = MatrixElement(
                    amptype, psp, tree=Process._tree_buf.value,
                    loop=LoopME(*Process._loop_buf),
                    iop=IOperator(*Process._ir1_buf),
                    loop2=Loop2ME(*Process._loop2_buf),
                    acc=Process._acc_buf.value)
        elif amptype == TREE:
            evaluate_tree_c(self.id, psp, Process._tree_buf)
            me = MatrixElement(amptype, psp, tree=Process._tree_buf.value)
        else:
            raise OpenLoopsError(
                'Process() amptype {} not implemented yet'.format(amptype))
        return me


class ProcessInfo(object):
    """Process info file data."""

    def __init__(self, libname):
        """Create ProcessInfo from info file for process 'libname'."""
        self.libname = libname
        self.channels = []
        self.mappings = {}
        self.condmappings = {}
        info_files = [fi for fi in os.listdir(proclib_dir)
                      if fi.endswith('.info')]
        info_files = [fi for fi in info_files
                      if fi.rsplit('_',1)[0] == 'libopenloops_' + libname]
        if not info_files:
            raise ProcessInfoError(
                'Process library {} info file not found.'.format(libname))
        self.loops_spec = info_files[0].rsplit('_',1)[1][:-5]
        with open(os.path.join(proclib_dir, info_files[0])) as fh:
            info = fh.readlines()
        for inf in info:
            if inf.startswith('options'): continue
            inf = inf.split()
            if inf[2].startswith('map='):
                # mapping[from] = to
                self.mappings[inf[1]] = inf[2].split('=')[1]
            elif inf[2].startswith('condmap='):
                # condmappings[from] = (to, [(param, val), ...])
                self.condmappings[inf[1]] = (
                    inf[2].split('=')[1], [cd.split('=') for cd in inf[3:]])
            else:
                # channels: (channel, number, {param: val), ...})
                self.channels.append(
                    (inf[1], int(inf[2]),
                     dict([cond.split('=') for cond in inf[3:]])))
        # sort by (loops_spec, channel, number)
        self.channels.sort(key=lambda ch:
                           (ch[2].get('Type', self.loops_spec), ch[0], ch[1]))

    def iterchannels(self, select=None):
        """Return a generator which provides channel identifiers of the form
        (<libname>_<channel>_<number>, loops_spec), taking into account
        conditional mappings with currently initialised parameters. I.e. only
        channels are returned which are not mapped to another channel.
        """
        # initialisation is needed to make get_parameter_double()
        # return the correct values
        parameters_flush()
        for ch in self.channels:
            cmap = self.condmappings.get(ch[0], None)
            do_map = False
            if cmap:
                # NOTE: only works with conditions
                #       for double precision parameters.
                do_map = all(get_parameter_double(param) == float(val)
                             for param, val in cmap[1])
            if not do_map:
                yield ('{}_{}_{}'.format(self.libname, ch[0], ch[1]),
                       ch[2].get('Type', self.loops_spec))

    def register(self, select=None):
        """Register channels which are returned by iterchannels()
        and return a list of Process() objects."""
        processes = []
        for ch in self.iterchannels(select):
            processes.append(Process(*ch))
        return processes


if __name__ == '__main__':
    set_parameter('stability_mode',11)
    proc = Process('1 -1 -> 22 21')
    me = proc.evaluate()
    print('pp =', me.psp)
    print(me)
