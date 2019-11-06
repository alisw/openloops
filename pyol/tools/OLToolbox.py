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


import os
import collections
import hashlib
import time
import sys

try:
    strtype = basestring
except NameError:
    strtype = str

timeformat = '%Y-%m-%d-%H-%M-%S'

def remove_duplicates(ls):
    seen = set()
    # seen.add() returns None
    return [el for el in ls if not (el in seen or seen.add(el))]

# =================================================== #
# functions to manage lists and dictionaries in files #
# =================================================== #


def strip_comments(ls):
    """Strip whitespace and comments starting with #
    in a list of strings and remove empty elements."""
    strippedls = [li.split('#')[0].strip() for li in ls]
    strippedls = [li for li in strippedls if li != '']
    return strippedls



def import_list(filename, lines=None, fatal=True,
                error_message='ERROR in import_list: file %s does not exist or is not readable.'):
    """Import a file as a list of strings stripped of whitespace, comments,
    and empty elements. filename can be a file name or a file object."""
    # OLD PYTHON
    decode = False
    if isinstance(filename, strtype):
        if not filename.startswith('http:'):
            try:
                fh = open(filename, 'r')
            except IOError:
                if fatal:
                    if '%s' in error_message:
                        print(error_message % (filename,))
                    else:
                        print(error_message)
                    raise
                else:
                    return None
        else:
            if sys.version_info < (3,0,0):
                from urllib2 import urlopen, URLError
            else:
                from urllib.request import urlopen
                from urllib.error import URLError
                decode = True
            try:
                fh = urlopen(filename)
            except URLError:
                if fatal:
                    if '%s' in error_message:
                        print(error_message % (filename,))
                    else:
                        print(error_message)
                    raise
                else:
                    return None
    else:
        fh = filename
    if lines:
        ls = [fh.readline() for n in range(lines)]
    else:
        ls = fh.readlines()
    fh.close()
    if decode:
        ls = [li.decode() for li in ls]
    return strip_comments(ls)



def split_to_pair(rec):
    """Split a string into a list with two elements at the first whitespace.

    If there is no whitespace, the second element will be an empty string."""
    pair = rec.split(None, 1)
    if len(pair) == 1:
        pair.append('')
    return pair


def split_to_dictionary(ls):
    """Split a list of strings into a dictionary at the first whitespace."""
    return dict([split_to_pair(rec) for rec in ls])



def import_dictionary(filename, fatal=True,
                      error_message='ERROR in import_dictionary: file %s does not exist or is not readable.'):
    """Read a dictionary from a file.

    Read the file with import_list() and create
    a dictionary with split_to_dictionary()."""
    dic = import_list(filename, fatal=fatal, error_message=error_message)
    if dic is None:
        return None
    else:
        return split_to_dictionary(dic)



def export_list(filename, ls):
    """Write a list of strings to a file, separated by newline."""
    #with open(filename, 'w'):
        #for el in ls:
            #fh.write(el)
    # OLD PYTHON
    try:
        fh = open(filename, 'w')
    except IOError:
        print('export_list: cannot open file', filename, 'for writing.')
        raise
    for el in ls:
        fh.write(el)
        fh.write('\n')
    fh.close()



def export_dictionary(filename, dic, form='%s %s'):
    """Export a dictionary to a file."""
    #ls = [form.format(key, val) for key, val in dic.items()]
    # OLD PYTHON
    ls = [form % (key, val) for key, val in dic.items()]
    export_list(filename, ls)


# ============ #
# git revision #
# ============ #

def get_git_revision(mandatory=False):
    """Get the git revision number from `git log -1`
    in the current working directory."""
    import subprocess
    git_exitcode = 1
    try:
        git_proc = subprocess.Popen(
            ['git', 'rev-parse', '--short', 'HEAD'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        git_out, git_err = git_proc.communicate()
        git_exitcode = git_proc.returncode
    except OSError:
        pass
    revision = 'none'
    if not git_exitcode:
        revision = git_out.decode().strip()
    if mandatory and (revision == 'none' or git_exitcode != 0):
        raise OSError(git_exitcode,
                      '`git info` failed. ' + git_err.decode().strip())
    return revision


# ============================ #
# Process library source files #
# ============================ #

def get_subprocess_src(loops, sub_process, processlib_src_dir, config,
                       nvirtualfiles=0, override_loops=False):
    """Return lists of double precision, multi precision and info files
    which belong to a subprocess. Used to determine which files must be
    compiled and which files are generated. Files to compile must be obtained
    with override_loops=True. Files to generate must be obtained with
    override_loops=False and the correct loops argument (taken from the
    process definition; the info file does not exist yet)."""
    info_files = [os.path.join(processlib_src_dir,
                               'info_' + sub_process + '.txt')]
    # Read info file for subprocess and
    # * check for 'Type' option to override loops specification,
    # * check for 'OLMode' option yo set olmode.
    olmode = 0
    if override_loops:
        try:
            fh = open(info_files[0], 'r')
        except IOError:
            print("Error reading process info file", info_files[0])
            raise
        info = fh.read()
        fh.close()
        for opt in info.split():
            if opt.startswith('Type='):
                loops = opt[5:]
            elif opt.startswith('OLMode='):
                olmode = int(opt[7:])

    if loops == 't':
        dp_src = ['born_generic_' + sub_process + '.F90']
        mp_src = ['born_' + sub_process + '.F90']
    elif 'l' in loops:
        if nvirtualfiles == 0:
            # determine the number of virtual files
            nvirtualfiles = int(import_list(os.path.join(
                  processlib_src_dir, 'virtual_' + sub_process + '.F90'),
                lines = 1)[0][1:])
        dp_src = ['loop_generic_' + sub_process + '.F90']
        mp_src = (['loop_' + sub_process + '.F90',
                   'virtual_' + sub_process + '.F90',
                   'tensorsum_' + sub_process + '.F90',
                   'checks_' + sub_process + '.F90']
                + ['virtual_' + str(n+1) + '_' + sub_process + '.F90'
                   for n in range(nvirtualfiles)])
        if 't' in loops:
            dp_src.append('born_generic4loop_' + sub_process + '.F90')
            mp_src.append('born4loop_' + sub_process + '.F90')
        if 'p' in loops:
            mp_src.append('pseudotree_' + sub_process + '.F90')
        # Decide if the loop process code is duplicated
        # to compiled it in quad precision.
        enable_process_qp = True
        if 's' in loops and not 't' in loops:
            # Loop-squared: qp only for checks
            enable_process_qp = config['process_qp_checks']
        else:
            if olmode <= 1:
                # OL1-type NLO or helicity optimised NLO: qp as rescue mode
                enable_process_qp = config['process_qp_rescue']
            else:
                # On-the-fly reduction: qp only for checks
                enable_process_qp = config['process_qp_checks']
        if not enable_process_qp:
            dp_src.extend(mp_src)
            mp_src = []

    dp_src = [os.path.join(processlib_src_dir, srcfile) for srcfile in dp_src]
    mp_src = [os.path.join(processlib_src_dir, srcfile) for srcfile in mp_src]

    return dp_src, mp_src, info_files



def get_processlib_src(loops, processlib, config):
    """Return lists of double precision, multi precision and info files
    which belong to a process library. Used to determine which files must be
    compiled. Not used to determine which files are generated."""
    processlib_src_dir = os.path.join(config['process_src_dir'], processlib)
    subprocesses = import_list(os.path.join(
        processlib_src_dir, 'process_definition', 'subprocesses.list'))
    subprocesses_extra = import_list(os.path.join(
        processlib_src_dir, 'process_definition', 'subprocesses_extra.list'),
        fatal=False)
    if subprocesses_extra and config['compile_extra']:
        subprocesses.extend(subprocesses_extra)
    dp_src = [os.path.join(processlib_src_dir,
                           'version_' + processlib + '.F90')]
    mp_src = []
    info_files = []

    for sub_process in subprocesses:
        dp_add, mp_add, info_add = get_subprocess_src(
            loops, sub_process, processlib_src_dir, config, override_loops=True)
        dp_src.extend(dp_add)
        mp_src.extend(mp_add)
        info_files.extend(info_add)

    return dp_src, mp_src, info_files



# ========================================= #
# Version database for a process repository #
# ========================================= #


class ProcessDB:
    """Process version database.

    Contains a dictionary with elements
    'process name' (string): 'version number' (integer)
    and optionally a file name from which the data was imported
    or to which the data was exported."""
    no_description = '"?|?|?|no description available"'

    def __init__(self, db=None, processes={}):
        """Create a process version database.

        Arguments:
        - db (optional) - the name of a process version database file to import.
        - processes (optional) - a dictionary
          {process: (date, hash, description), ...}.
        First the database is imported, then it is updated with 'processes'."""
        if db:
            self.import_db(db)
        else:
            self.content = {}
            self.db_file = None
        self.updated = False
        self.update(processes)

    def update(self, processes):
        """Update versions of existing processes or add new processes.
        Set the process hash and upload date in the database.

        Arguments:
        - processes: see __init__."""
        if processes:
            self.content.update(processes)
            self.updated = True

    def import_db(self, db):
        """Import process version database from a file."""
        self.db_file = db
        data = import_dictionary(db)
        # OLD PYTHON
        self.content = dict([(proc, (lambda ls: ls + [self.no_description]
                                     if len(ls) == 2 else ls
                                     )(hash_date_descr.split(None,2)))
                             for proc, hash_date_descr in data.items()])

    def remove(self, processes):
        """Remove processes from a process version database."""
        if not isinstance(processes, (list, tuple, set)):
            processes = [processes]
        old_size = len(self.content)
        for proc in processes:
            self.content.pop(proc, None)
        if len(self.content) != old_size:
            self.updated = True

    def export_db(self, db=None):
        """Export process version database to a file."""
        if self.updated or (db is not None and db != self.db_file):
            if db:
                self.db_file = db
            tmp_file = self.db_file + '.~' + str(os.getpid())
            tmp_content = dict([(key, '  '.join(val))
                                for key, val in self.content.items()])
            export_dictionary(tmp_file, tmp_content, '%-30s %s')
            os.rename(tmp_file, self.db_file)

    def __getitem__(self, key):
        return self.content[key]

    def get(self, key, default=None):
        return self.content.get(key, default)



class ChannelDB:
    """Database for partonic channels with process info."""

    def __init__(self, db=None, channels={}):
        """Create a channel database.

        Arguments:
        - db (optional) - the name of a channel database file to import.
        - channels (optional) - library channel info to add to the database,
          dictionary {libname: [channel_info_1, ...], ...}
        First the database is imported, then it is updated with 'channels'."""
        if db:
            self.import_db(db)
        else:
            self.content = {}
            self.db_file = None
        self.updated = False
        self.update(channels)

    def update(self, channels):
        """Update the channel database: delete old channels
        contained in libraries whose names are in 'channels'.

        Arguments:
        - channels: see __init__."""
        if channels:
            self.content.update(channels)
            self.updated = True

    def import_db(self, db):
        """Import channel database from a file."""
        self.db_file = db
        data = [ch.split() for ch in import_list(db)[1:]]
        self.content = collections.defaultdict(list)
        for iline in data:
            self.content[iline[0]].append(iline)
        self.content = dict(self.content)

    def remove(self, libnames):
        """Remove channels from the database
        which belong to a library in 'libnames'."""
        if not isinstance(libnames, (list, tuple, set)):
            libnames = [libnames]
        old_size = len(self.content)
        for proc in libnames:
            self.content.pop(proc, None)
        if len(self.content) != old_size:
            self.updated = True

    def export_db(self, db=None):
        """Export channel database to a file."""
        if self.updated or (db is not None and db != self.db_file):
            if db:
                self.db_file = db
            tmp_file = self.db_file + '.~' + str(os.getpid())
            data = []
            for proc in sorted(self.content.keys()):
                data.extend([' '.join(ch) for ch in self.content[proc]])
            channels_hash = hashlib.md5()
            for iline in data:
                channels_hash.update(iline.encode())
            data.insert(0, channels_hash.hexdigest() + '  ' +
                           time.strftime(timeformat))
            export_list(tmp_file, data)
            os.rename(tmp_file, self.db_file)


def repo_name(repo):
    # (assumes that there are no public repositories named .+_.{16})
    if len(repo) > 16 and repo[-17] == '_' and '_' not in repo[-16:]:
        return repo[:-17]
    else:
        return repo
