#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

# This module requires Python 2.7, but this must be tested before importing
# the module.

import sys
import os
import shlex
import math
import atexit
import pickle
from ctypes import CDLL, byref, c_int, c_double, c_char_p

import OLBaseConfig

if sys.platform.startswith("darwin"):
    shared_lib_ext = ".dylib"
else:
    shared_lib_ext = ".so"

pyol = CDLL("libopenloops" + shared_lib_ext)

pyol_setparameter = pyol.ol_setparameter_string
pyol_rambo = pyol.ol_rambo
pyol_parameters_flush = pyol.ol_parameters_flush
pyol_parameters_write = pyol.ol_printparameter
pyol_finish = pyol.ol_finish

pyol_n_particles = "ol_n_external_%s"
pyol_get_external_masses = "ol_get_masses_%s"
pyol_set_permutation = "ol_set_permutation_%s"
pyol_amp2tree = "ol_amp2tree_%s"
pyol_vamp2 = "ol_vamp2_%s"
pyol_ctamp2 = "ol_ctamp2_%s"
pyol_ptamp2 = "ol_ptamp2_%s"


if sys.version_info[:2] < (2,7):
    print "This module requires Python 2.7."
    sys.exit(1)



# Openloops interface version.
# In future version this can be used to import ProcessTestData in a deprecated format.
version = (0,1)
test_data_folder = "test_data"
# center-of-mass energy if not specified otherwise
default_energy = 1000
# number of phase space points in test data if not specified otherwise
default_sample_size = 10
# various accuracy thresholds in digits
phasespace_accuracy_warn = 9
phasespace_accuracy_fatal = 5
tree_accuracy_warn = 9
tree_accuracy_fatal = 5
matrix_element_agreement = 5
# precision of 'float' in digits
max_precision = sys.float_info.mant_dig * math.log(2,10)


atexit.register(pyol_finish)


class OpenLoopsError(Exception):
    pass



def float_list_string(ls):
    return "[ " + ", ".join(["{:>23.15e}".format(val) for val in ls]) + " ]"



def common_digits(x, y):
    if x == y:
        return max_precision
    else:
        return -math.log(2*abs(x-y)/(abs(x)+abs(y)), 10)



def valid_permutation(arg, n=None):
    if isinstance(arg, str):
        try:
            perm = eval(arg)
        except SyntaxError:
            raise ValueError("Invalid permutation: " + repr(arg))
    else:
        perm = arg
    if not (isinstance(perm, list) and
            all([lambda nn: isinstance(nn, int) for nn in perm]) and
            set(perm) == set(range(1, len(perm)+1))):
        raise ValueError("Invalid permutation: " + repr(arg))
    if n is not None and max(perm) != n:
        raise ValueError("Permutation " + repr(arg) + " does not have length " + repr(n))
    return perm



def to_subprocess_id(sp):
    if isinstance(sp, SubprocessID):
        return sp
    else:
        return SubprocessID(sp)



def parse_subprocess_argument(arg, process=None):
    """
    Convert a string into a list of SubprocessID instances.

    Example: "sub[2,3,1],sub_2,..." --> [SubprocessID("sub[2,3,1]"), SubprocessID("sub_2"), ...]
    """
    # quote permutation lists: sub[1,2,3] --> sub"[1,2,3]"
    if isinstance(arg, str):
        sl = arg.replace("[", "\"[").replace("]", "]\"")
        sl = shlex.shlex(sl, posix=True)
        sl.whitespace = ","
        sl.whitespace_split = True
        sl.wordchars += "[]"
        sp_list = [SubprocessID(sp) for sp in sl]
    elif isinstance(arg, (list, tuple, set)):
        sp_list = [to_subprocess_id(sp) for sp in arg]
    else:
        raise ValueError("Invalid subprocess argument: " + repr(arg))
    if process is not None:
        sp_list = [sp.validate(process) for sp in sp_list]
    return sp_list



class Parameters:
    """
    Parameters class to store, and commit parameters via pyol_setparameter().

    No default values are provided, parameters must be reset explicitly if needed
    when they were initialised before, e.g. by a different Parameters instance.
    """
    def __init__(self, commit = True, **parameters):
        self.parameters = {key: str(val) for key, val in parameters.items()}
        if commit:
            for key, val in self.parameters.items():
                pyol_setparameter(c_char_p(key), c_char_p(val))
            pyol_parameters_flush()

    set = __init__

    @staticmethod
    def write():
        """Call the Fortran parameters_write() subroutines."""
        pyol_parameters_write()

    def __repr__(self):
        return "Parameters(" + ", ".join([" = ".join(pair) for pair in self.parameters.items()]) + ")"

    def __str__(self):
        return "\n".join([" = ".join(pair) for pair in self._parameters.items()])



class PhaseSpacePoint:
    """
    A phase space point consisting of masses, energy, and momenta.

    Arguments:
    masses -- list of masses or the number of massless particles
    energy -- center of mass energy; has default default_energy
    momenta -- phase space point as c_double array;
               if momenta is present, masses and energy must not be present
    """
    @staticmethod
    def mass_from_momentum(momentum):
        mass = momentum[0]**2 - momentum[1]**2 - momentum[2]**2 - momentum[3]**2
        if abs(mass) < phasespace_accuracy_warn:
            mass = 0.
        # Commented so that inconsistencies are handled within Fortran
        #elif mass < 0:
            #raise OpenLoopsError("mass < 0 in phase space piont")
        return math.sqrt(mass)

    def __init__(self, masses=None, energy=None, momenta=None):
        if momenta is None:
            # Generate a phase space point with RAMBO
            if energy is None:
                energy = default_energy
            if masses is None:
                raise OpenLoopsError("masses needed")
            elif isinstance(masses, int) and masses >= 4:
                masses = masses * [0.]
            self.masses = list(masses) # copy
            self.energy = energy

            n_particles = len(masses)
            c_masses_type = n_particles * c_double
            c_psp_type = n_particles * (4*c_double)
            p_rambo = c_psp_type()

            pyol_rambo(byref(c_double(energy)),
                       byref(c_masses_type(*masses)),
                       byref(c_int(n_particles)),
                       byref(p_rambo))

            self.momenta = p_rambo

        else:
            # A phase space point is given
            if isinstance(momenta, PhaseSpacePoint):
                self.masses = momenta.masses
                self.energy = momenta.energy
                self.momenta = momenta.momenta
            else:
                self.masses = map(PhaseSpacePoint.mass_from_momentum, momenta)
                self.energy = sum([mom[0] for mom in momenta])/2
                if isinstance(momenta, (list, tuple)):
                    # convert to a c_double array
                    self.momenta = (len(momenta) * (4*c_double))(*map(tuple, momenta))
                else:
                    # should already be a c_double array
                    self.momenta = momenta
            consistency = 1000
            if masses is not None:
                consistency = min([common_digits(m1,m2) for m1, m2 in zip(*[self.masses, masses])])
            if energy is not None:
                consistency = min(consistency, common_digits(self.energy, energy))
            if consistency < phasespace_accuracy_warn:
                if consistency < phasespace_accuracy_fatal:
                    raise OpenLoopsError("inconsistent with given masses/energy, agreement: "
                      + str(consistency) + " digits")
                else:
                    print ("WARNING: phase space point is compatible with given masses/energy only to "
                      + str(consistency) + " digits.")

    def __repr__(self):
        return "PhaseSpacePoint(masses = {}, energy = {}, momenta =\n{})".format(
               self.masses, self.energy, map(list, self.momenta))

    def __str__(self):
        return str(map(list, self.momenta))



class Process:
    """
    Process class.
    Provides access to matrix element routines.
      self.permute[subprocess]([permutation]) = active permutation
      self.external[subprocess]() = number_of_particles, masses
      self.tree[subprocess](ps_point) = M2tree
      self.loop[subprocess](ps_point) = M2tree, M2loop0, M2loop1, M2loop2, IR0, IR1, IR2
      self.ct[subprocess](ps_point) = M2tree, M2ct
      self.pseudotree[subprocess](ps_point) = M2tree, M2pt, M2loop
    """
    _proclib_path = "proclib"
    _m2l0 = c_double()
    _m2ct = c_double()
    _m2pt = c_double()
    _m2l1 = (3*c_double)()
    _irl1 = (3*c_double)()
    _m2l2 = (5*c_double)()
    _irl2 = (5*c_double)()

    @staticmethod
    def _get_library_name(name, loops="auto"):
        """
        Check for existing libraries for a process and return a dictionary with their names.

        Arguments:
          name -- the process name
          loop -- the loops specification (defaults to "auto")
        If loops="auto", the call only succeeds if the library is unique.
        The result conatins the keys "loops", "process_module" (file name), and "info_file" (file name).
        """
        process_library = "libopenloops_" + name + "_%s" + shared_lib_ext
        info_file = "libopenloops_" + name + "_%s.info"

        if loops == "auto":
            available_loops = [lps for lps in OLBaseConfig.loops_specifications
                               if lps != "auto" and os.path.isfile(os.path.join(Process._proclib_path,
                                    process_library % (lps,)))]
            if len(available_loops) == 0:
                raise OpenLoopsError("process library not found " + repr(name))
            elif len(available_loops) > 1:
                raise OpenLoopsError("process library not unique " + repr(available_loops))
            loops = available_loops[0]
        else:
            if loops not in OLBaseConfig.loops_specifications:
                raise OpenLoopsError("unknown loops specification " + repr(loops))
            if not os.path.isfile(os.path.join(Process._proclib_path, process_library % (loops,))):
                raise OpenLoopsError("process library not found: " + repr(name) + "_" + repr(loops))
        if not os.path.isfile(os.path.join(Process._proclib_path, info_file % (loops,))):
            raise OpenLoopsError("info file not found: " + repr(name) + "_" + repr(loops))
        return {"loops": loops,
                "process_library": process_library % (loops,),
                "info_file": info_file % (loops,)}

    @staticmethod
    def _parse_infoline(infoline):
        """Convert a string "KEY1=VAL1 KEY2=VAL2 ..." to a dictionary {KEY1: VAL1, KEY2: VAL2, ...}."""
        il = infoline.split()
        return [il[1], {il[2]: dict([[opt[0], opt[1]] for opt in [opt.split("=") for opt in il[3:]]])}]
        #return (il[1] + "_" + il[2], dict([[opt[0], eval(opt[1])] for opt in [opt.split("=") for opt in il[3:]]]))

    @staticmethod
    def _get_info(name=None, loops="auto", info_file=None):
        """
        Read an info file and return a dictionary with the available sub-processes.

        Arguments:
          name -- the process name
          loops -- the loops specification (defaults to "auto")
          info_file -- name of the info file
        Either name or info_file must be given.
        This is a static method so that it can be called without creating a Process instance.
        """
        if info_file is None:
            info_file = Process._get_library_name(name, loops)["info_file"]
        if info_file == os.path.basename(info_file):
            info_file = os.path.join(Process._proclib_path, info_file)
        # Fill a dictionary with
        # {SubProcessA_1: {VarA11: ValA11, ...}, SubProcessA_2: {...}, SubProcessB_1: ..., ...}
        available_subprocesses = {}
        with open(info_file) as fh:
            for infoline in fh.readlines():
                if " map=" in infoline or " condmap=" in infoline or infoline.startswith("options "):
                    continue
                subprocinfo = Process._parse_infoline(infoline)
                subprocnum = subprocinfo[1].keys()[0] # there is only one element here
                available_subprocesses[subprocinfo[0] + '_' + str(subprocnum)] = subprocinfo[1][subprocnum]
        return dict(available_subprocesses)
        # Fill a dictionary with
        # {SubProcessA_1: {VarA11: ValA11, ...}, SubProcessA_2: {...}, SubProcessB_1: {...}, ...}
        #available_subprocesses = {}
        #with open(info_file) as fh:
            #for infoline in fh.readlines():
                #subprocinfo = Process._parse_infoline(infoline)
                #available_subprocesses[subprocinfo[0]] = subprocinfo[1]
        #return available_subprocesses


    def __init__(self, name, loops="auto"):
        library_name = Process._get_library_name(name, loops)
        loops = library_name["loops"]
        process_library = library_name["process_library"]

        self.name = name        # process name
        self.loops = loops      # loops specification
        self.info = Process._get_info(name, loops)
        self.subprocesses = sorted(self.info.keys())
        self.external = {}      # dictionary of functions
        self.nexternal = {}     # number of external particles
        self.permute = {}       # dictionary of functions
        self.rambo = {}         # dictionary of functions
        self.tree = {}          # dictionary of functions
        self.loop = {}          # dictionary of functions
        self.ct = {}            # dictionary of functions
        self.pseudotree = {}    # dictionary of functions
        self.matrix_element = {} # dictionary pointing to the MatrixElement constructor

        self._library = CDLL(process_library)
        self._fortran_external = {}
        self._fortran_permutation = {}
        self._active_permutation = {}
        self._fortran_tree = {}
        self._fortran_loop = {}
        self._fortran_ct = {}
        self._fortran_pseudotree = {}

        self.contains_tree = {}
        self.contains_loop = {}
        self.contains_pseudotree = {}
        self.contains_loopsquared = {}

        if "s" in loops:
            self.loop_type = "loopsquared"
        elif "l" in loops:
            self.loop_type = "loop"
        else:
            self.loop_type = "tree"

        for subprocess in self.subprocesses:
            process_name = (self.name + "_" + subprocess).lower()

            force_loops = self.info[subprocess].get('Type', None)
            if not force_loops:
                force_loops = loops
            self.contains_tree[subprocess] = "t" in force_loops
            self.contains_loop[subprocess] = "l" in force_loops
            self.contains_pseudotree[subprocess] = "p" in force_loops
            self.contains_loopsquared[subprocess] = "s" in force_loops

            n_ext = c_int()
            getattr(self._library, pyol_n_particles % (process_name,))(byref(n_ext))
            self.nexternal[subprocess] = n_ext.value
            self._active_permutation[subprocess] = range(1, self.nexternal[subprocess] + 1)

            self._fortran_external[subprocess] = getattr(self._library, pyol_get_external_masses % (process_name,))

            def external_tmp(subprocess=subprocess):
                """Get the masses of the external particles in the active permutation."""
                c_masses = (self.nexternal[subprocess]*c_double)()
                self._fortran_external[subprocess](byref(c_masses))
                return list(c_masses)

            self.external[subprocess] = external_tmp
            del external_tmp

            self._fortran_permutation[subprocess] = getattr(self._library, pyol_set_permutation % (process_name,))

            def permutation_tmp(permutation=None, subprocess=subprocess):
                """
                Set the permutation for the sub-process if it changed and return it.
                If no argument is given, the active permutation is returned.
                """
                if permutation is not None and permutation != self._active_permutation[subprocess]:
                    c_perm_type = self.nexternal[subprocess] * c_int
                    self._fortran_permutation[subprocess](byref(c_perm_type(*permutation)))
                    self._active_permutation[subprocess] = list(permutation)
                return self._active_permutation[subprocess]

            self.permute[subprocess] = permutation_tmp
            del permutation_tmp

            def rambo_tmp(energy_or_momenta=None, subprocess=subprocess):
                """
                Call the PhaseSpacePoint constructor for the sub-process.

                If no argument is given, None will be used (resulting in default_energy).
                If the argument is a phase space point (as a list/array, or as a PhaseSpacePoint instance),
                return it as PhaseSpacePoint instance.
                """
                if energy_or_momenta is None or isinstance(energy_or_momenta, (int, float)):
                    return PhaseSpacePoint(self.external[subprocess](), energy_or_momenta)
                else:
                    return PhaseSpacePoint(momenta=energy_or_momenta)

            self.rambo[subprocess] = rambo_tmp
            del rambo_tmp

            if self.contains_tree[subprocess]:
                self._fortran_tree[subprocess] = getattr(self._library, pyol_amp2tree % (process_name,))
                def tree_tmp(ps_pt=None, permutation=None, subprocess=subprocess):
                    """
                    Calculate a tree matrix element. Return M2tree.

                    Arguments:
                      ps_pt -- either a phase space point (list/array/PhaseSpacePoint),
                               or the energy, or None (defaults to None, resulting in default_energy)
                      permutation -- the permutation of the external particles (optional)
                    """
                    self.permute[subprocess](permutation)
                    ps_pt = self.rambo[subprocess](ps_pt)

                    self._fortran_tree[subprocess](byref(ps_pt.momenta), byref(self._m2l0))
                    # return a tuple with one element
                    return (self._m2l0.value,)

                self.tree[subprocess] = tree_tmp
                del tree_tmp

            if self.contains_loop[subprocess]:
                self._fortran_loop[subprocess] = getattr(self._library, pyol_vamp2 % (process_name,))
                def loop_tmp(ps_pt=None, permutation=None, subprocess=subprocess):
                    """
                    Calculate a loop matrix element. Return (M2L0, M2L1(0:2), IRL1(0:2), M2L2(0:4), IRL2(0:4)).

                    M2L1(0)..M2L1(2) are the finite, 1/ep and 1/ep^2 parts of the loop matrix element.
                    IRL1(0)..IRL1(2) are the finite, 1/ep and 1/ep^2 parts of the i-operator.
                    M2L2(0)..M2L2(4) and IRL2_0..IRL2_4 are the loop^2 matrix element and IR contribution
                    as a Laurant series.

                    Arguments:
                      ps_pt -- either a phase space point (list/array/PhaseSpacePoint),
                               or the energy, or None (defaults to None, resulting in default_energy)
                      permutation -- the permutation of the external particles (optional)
                    """
                    self.permute[subprocess](permutation)
                    ps_pt = self.rambo[subprocess](ps_pt)

                    self._fortran_loop[subprocess](byref(ps_pt.momenta), byref(self._m2l0),
                                                   byref(self._m2l1), byref(self._irl1),
                                                   byref(self._m2l2), byref(self._irl2))

                    return (self._m2l0.value, list(self._m2l1), list(self._irl1),
                                              list(self._m2l2), list(self._irl2))

                self.loop[subprocess] = loop_tmp
                del loop_tmp

                self._fortran_ct[subprocess] = getattr(self._library, pyol_ctamp2 % (process_name,))
                def ct_tmp(ps_pt=None, permutation=None, subprocess=subprocess):
                    """
                    Calculate a counter-term matrix element. Return (M2tree, M2ct).

                    Arguments:
                      ps_pt -- either a phase space point (list/array/PhaseSpacePoint),
                               or the energy, or None (defaults to None, resulting in default_energy)
                      permutation -- the permutation of the external particles (optional)
                    """
                    self.permute[subprocess](permutation)
                    ps_pt = self.rambo[subprocess](ps_pt)

                    self._fortran_ct[subprocess](byref(ps_pt.momenta), byref(self._m2l0), byref(self._m2ct))

                    return (self._m2l0.value, self._m2ct.value)

                self.ct[subprocess] = ct_tmp
                del ct_tmp

            if self.contains_pseudotree[subprocess]:
                self._fortran_pseudotree[subprocess] = getattr(self._library, pyol_ptamp2 % (process_name,))
                def pseudotree_tmp(ps_pt=None, permutation=None, subprocess=subprocess):
                    """
                    Calculate a pseudo-tree matrix element. Return (M2tree, M2pt, M2loop).

                    Arguments:
                      ps_pt -- either a phase space point (list/array/PhaseSpacePoint),
                               or the energy, or None (defaults to None, resulting in default_energy)
                      permutation -- the permutation of the external particles (optional)
                    """
                    self.permute[subprocess](permutation)
                    ps_pt = self.rambo[subprocess](ps_pt)

                    self._fortran_pseudotree[subprocess](byref(ps_pt.momenta), byref(self._m2l0),
                                                         byref(self._m2pt), byref(self._m2ct)) # loop --> ct

                    return (self._m2l0.value, self._m2pt.value, self._m2ct.value)

                self.pseudotree[subprocess] = pseudotree_tmp
                del pseudotree_tmp

            def matrix_element_tmp(ps_point=None, permutation=None, subprocess=subprocess):
                """
                Call the MatrixElement constructor.
                """
                return MatrixElement(self, SubprocessID(subprocess, permutation), ps_point)
            self.matrix_element[subprocess] = matrix_element_tmp
            del matrix_element_tmp

    def __repr__(self):
        return "Process(name={}, loops={})".format(self.name, self.loops)

    def __str__(self):
        return "process name:  {}\nloops:         {}\nsub-processes: {}\n".format(
               self.name, self.loops, ", ".join(self.subprocesses))



class SubprocessID:

    def __init__(self, subprocess, permutation=None):
        if not isinstance(subprocess, str):
            raise ValueError("Invalid sub-process identifier: " + repr(subprocess))
        subprocess = subprocess.split("[", 1)
        permutation_from_subprocess = None
        if len(subprocess) == 2:
            permutation_from_subprocess = valid_permutation("[" + subprocess[1])
        subprocess = subprocess[0]

        if not subprocess.split("_")[-1].isdigit():
            subprocess = subprocess + "_1"
        if permutation is None:
            permutation = permutation_from_subprocess
        else:
            valid_permutation(permutation)
            if permutation_from_subprocess is not None:
                raise ValueError("Permutation cannot be given twice for a sub-process.")

        self.name = subprocess
        self.permutation = permutation

    def validate(self, process, default_permutation=None):
        # Check if self.name is in process.subprocesses
        # and if set(self.permutation) == set(range(1, n_particles+1)).
        # If self.permutation == None, replace it by default_permutation,
        # if present, otherwise range(1, n_particles+1).
        if not isinstance(process, Process):
            raise ValueError("process must be a process instance.")
        if not self.name in process.subprocesses:
            raise OpenLoopsError("Process " + repr(process.name)
              + " does not contain a sub-processes " + repr(self.name))
        if self.permutation is None:
            if default_permutation is None:
                self.permutation = range(1, process.nexternal[self.name] + 1)
            else:
                self.permutation = valid_permutation(default_permutation)
        else:
            if set(self.permutation) != set(range(1, process.nexternal[self.name] + 1)):
                raise ValueError("Invalid permutation: " + str(self.permutation) + " for sub-process " + self.name)
        return self

    def __str__(self):
        if self.permutation is None:
            perm = ""
        else:
            perm = str(self.permutation).replace(" ", "")
        return self.name + perm

    def __repr__(self):
        return "SubprocessID(" + repr(self.name) + ", " + str(self.permutation).replace(" ", "") + ")"



class MatrixElement:
    """
    Store a phase space point and the matrix elements at this point for a partonic process.

    Members:
      ps_point -- a PhaseSpacePoint instance
      value -- the matrix elements as returned by Fortran
               (float for tree process, 7-tuple for loop)

    Does NOT contain any further information like the process, parameters, permutation, ...
    """
    def __init__(self, process=None, subprocess=None, ps_point=None,
                 tree=True, loop=True, value=None, verbose=0):
        """
        Calculate matrix a element.

        Arguments:
          process -- a Process instance
          subprocess -- a (permuted) sub-process as SubprocessID() or as a string;
                        if omitted or None and there is only one sub-process, use this one
          ps_point -- if omitted or None, use default_energy,
                      the energy if it is a number,
                      or a phase space point, either as an array, or as a PhaseSpacePoint instance
          tree -- calculate the tree matrix element, if a tree routine is available
          loop -- calculate the loop matrix element, if a loop routine is available
          value -- set the matrix element to value instead of calculating it;
                   only makes sense if ps_point is actually a phase space point.
        """
        if process is not None:
            if subprocess is None:
                if len(process.subprocesses) == 1:
                    subprocess = list(process.subprocesses)[0]
                else:
                    raise OpenLoopsError("sub-process not unique")
            if not isinstance(subprocess, SubprocessID):
                subprocess = SubprocessID(subprocess)
            subprocess.validate(process)
            permutation = subprocess.permutation
            subprocess_name = subprocess.name

            looptree = tree and not loop and not process.contains_tree[subprocess_name]
            tree = tree and process.contains_tree[subprocess_name]
            loop = loop and process.contains_loop[subprocess_name]

            process.permute[subprocess_name](permutation)
            self.ps_point = process.rambo[subprocess_name](ps_point)
            if loop:
                self.value = process.loop[subprocess_name](self.ps_point)
                if tree:
                    tree = process.tree[subprocess_name](self.ps_point)
                    consistency = common_digits(self.value[0], tree[0])
                    if consistency < tree_accuracy_warn:
                        if consistency < tree_accuracy_fatal:
                            raise OpenLoopsError("inconsistent tree matrix element from tree and loop routines, agreement: " + repr(consistency) + " digits.")
                        else:
                            print ("WARNING: tree matrix element from tree and loop routines only consistent to"
                            + str(consistency) + " digits.")
                if process.contains_loopsquared[subprocess_name]:
                    # flatten tuple
                    self.value = (self.value[0],) + tuple(self.value[1]) + tuple(self.value[2]) + tuple(self.value[3]) + tuple(self.value[4])
                else:
                    # flatten tuple and discard loop^2 information
                    self.value = (self.value[0],) + tuple(self.value[1]) + tuple(self.value[2])
            elif tree:
                self.value = process.tree[subprocess_name](self.ps_point)
            elif looptree:
                self.value = (process.loop[subprocess_name](self.ps_point)[0],)
            else:
                raise OpenLoopsError("tree/loop request is not satisfiable.")
        elif process is None and value is not None:
            self.ps_point = ps_point
            self.value = value
        else:
            OpenLoopsError("Process or value needed to initialise matrix element.")
        if verbose > 1:
            print self.ps_point
        if verbose > 0:
            print float_list_string(self.value)

        #if process is not None and value is not None:
            #value_is_loop = isinstance()

    def compare(self, other):
        if len(self.value) == len(other.value):
            agreement = [common_digits(*vals) for vals in zip(self.value, other.value)]
        else:
            raise OpenLoopsError("Cannot compare matrix elements of different types.")
        return min(agreement)

    def __repr__(self):
        return "MatrixElement(ps_point=\n{!r},\nvalue={})".format(self.ps_point, self.value)

    def __str__(self):
        return str(self.value)



class ProcessTestData:
    """
    Proecss test data: store a sample matrix elements for a process
    and the environment that was used to generate them.

    Members:
      version -- version information as a tuple of two integers.
      process_name -- the process name.
      loop_type -- "tree", "loop", or "loopsquared".
      parameters -- the Parameters instance used to calculate the matrix elements.
      subprocess -- a list SubprocessID instances.
      all_subprocesses -- True if all sub-processes of the process are contained.
      matrix_elements -- a list of MatrixElement instances
    Temporary members:
      check_successful: set by validate(), deleted by dump()/dumps().
    """
    def __init__(self, process, subprocess=None, phasespace=None, n=None,
                 parameters=None, verbose=0):
        """
        Create sample data for a process.

        subprocess -- a list of (permuted) sub-processes as a string,
                      or a list of strings (possibly as a string)
                      or SubprocessID instances.
        phasespace -- a list of phase space points (as accepted by PhaseSpacePoint),
                      or a list of energies, or the energy.
        n -- number of matrix elements per sub-process.
        parameters -- the Parameters instance to be used.
        verbose -- 0: print nothing, 1: print matrix elements, 2: also pront phase space points.
        """
        self.version = version
        self.matrix_elements = {}

        if isinstance(process, str):
            process = Process(process)

        if n is None:
            n = default_sample_size

        if process is None:
            pass
        else:
            self.process_name = process.name
            self.loop_type = process.loop_type
            if parameters is None:
                self.parameters = Parameters()
            else:
                self.parameters = parameters
            self.all_subprocesses = subprocess is None

            self.parameters.set()

            if self.all_subprocesses:
                subprocess = process.subprocesses

            self.subprocesses = parse_subprocess_argument(subprocess, process)

            if phasespace is None or isinstance(phasespace, (int, float)):
                phasespace = n * [phasespace]

            for subprocess in self.subprocesses:
                if verbose > 0:
                    print "sub-process:", str(subprocess)
                self.matrix_elements[subprocess.name] = [
                    MatrixElement(process, subprocess, ps_point, verbose=verbose)
                    for ps_point in phasespace]
                if verbose > 0:
                    print

    def dump(self, filename=None):
        """Write a pickled representation to the file 'filename' (default 'processname.pkl')"""
        try:
            del self.check_successful
        except AttributeError:
            pass
        if filename is None:
            filename = self.process_name + ".pkl"
        try:
            os.remove(filename)
        except OSError:
            pass
        dump_dir = os.path.dirname(filename)
        if not os.path.exists(dump_dir):
            os.mkdir(dump_dir)
        with open(filename, "w") as fh:
            pickle.dump(self, fh)

    def dumps(self):
        """Return a pickled representation as a string."""
        try:
            del self.check_successful
        except AttributeError:
            pass
        return pickle.dumps(self)

    @staticmethod
    def load(filename):
        """Load pickled ProcessTestData from a file."""
        with open(filename, "r") as fh:
            ptd = pickle.load(fh)
        return ptd

    def validate(self, new_process_name=None, verbose=0):
        """
        Validate ProcessTestData.

        Calculates matrix elements with the settings stored in 'self'
        and compares the results to those in 'self'.

        new_process_name: set if the name of the process changed.
        verbose: see __init__.
        """

        agrees = True

        if new_process_name is None:
            new_process_name = self.process_name

        process = Process(new_process_name)
        new_ptd = ProcessTestData(process=None)

        new_ptd.process_name = new_process_name
        new_ptd.loop_type = process.loop_type
        new_ptd.parameters = self.parameters
        new_ptd.all_subprocesses = self.all_subprocesses

        if new_ptd.loop_type != self.loop_type:
            raise OpenLoopsError("ProcessTestData.validate(): loops regression in " + repr(new_ptd.process_name))

        new_ptd.parameters.set()

        self_sp_names = set([sp.name for sp in self.subprocesses])

        if self.all_subprocesses:
            new_ptd.subprocesses = parse_subprocess_argument(process.subprocesses, process)
            missing_sps = self_sp_names - set(process.subprocesses)
            new_sps = set(process.subprocesses) - self_sp_names
            if missing_sps:
                print "WARNING: subprocesses were removed from", self.process_name + ":", missing_sps
            if new_sps:
                print "WARNING: subprocesses were added to", self.process_name + ":", new_sps
        else:
            missing_sps = self_sp_names - set(process.subprocesses)
            if missing_sps:
                print "WARNING: subprocesses were removed from", self.process_name + ":", missing_sps
            new_ptd.subprocesses = [sp for sp in self.subprocesses if sp.name not in missing_sps]

        for subprocess in new_ptd.subprocesses:
            if verbose > 0:
                print "sub-process:", str(subprocess)
            subprocess_name = subprocess.name
            new_ptd.matrix_elements[subprocess_name] = []
            me_count = 0
            for me in self.matrix_elements[subprocess_name]:
                new_me = MatrixElement(process, subprocess, me.ps_point, verbose=verbose)
                if me.compare(new_me) < matrix_element_agreement:
                    if verbose > 0:
                        print float_list_string(me.value), "(old)"
                    print "WARNING: {}_{}: Matrix element {} disagrees:\n{}\n{}".format(
                        new_ptd.process_name, subprocess, me_count+1, me, new_me)
                    agrees = False
                new_ptd.matrix_elements[subprocess_name].append(new_me)
                me_count = me_count + 1
            if verbose > 0:
                print

        new_ptd.check_successful = agrees

        return new_ptd

    def show(self, verbose=1):
        """
        Print ProcessTestData.

        verbose -- 0: only process, loop type, version
                   1: also matrix elements
                   2: also phase space points
                   3: also used parameters
        """
        print "Process:", self.process_name, "loop_type", self.loop_type, "version", self.version
        print
        if verbose > 2:
            print self.parameters
            print
        for subprocess in self.subprocesses:
            if verbose > 0:
                print "sub-process:", str(subprocess)
                for me in self.matrix_elements[subprocess.name]:
                    if verbose > 1:
                        print me.ps_point
                    print float_list_string(me.value)
                print

    #def __repr__(self):

    #def __str__(self):
