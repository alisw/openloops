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
import os
import sys
import shutil
import argparse
import time
import hashlib

sys.path.insert(0, os.path.abspath(os.path.join('pyol', 'config')))
sys.path.insert(0, os.path.abspath(os.path.join('pyol', 'tools')))

import OLBaseConfig
import OLToolbox

if sys.version_info < (3,0,0):
    from urllib2 import urlopen, URLError
else:
    from urllib.request import urlopen
    from urllib.error import URLError

commandline_options = [arg.split('=',1) for arg in sys.argv[1:] if ('=' in arg and not arg.startswith('-'))]
config = OLBaseConfig.get_config(commandline_options)

if config['print_python_version']:
    print('download_process.py uses Python', sys.version)

parser = argparse.ArgumentParser(
    description="""Download process code from a web server if a newer version
                   is available which is compatible with the installed
                   OpenLoops version.""")
parser.add_argument('processes', metavar='proc', nargs='*',
                    help="""process to download,
                            treated as a collection if it ends with '/'.""")
parser.add_argument('-i', '--ignore', action='store_true',
                    help='ignore non-existing processes and collections')
parser.add_argument('-f', '--force', action='store_true',
                    help='force download')
args = parser.parse_args(
    [arg for arg in sys.argv[1:] if (arg.startswith('-') or '=' not in arg)])
procs = args.processes

print('\n>>> OpenLoops Process Downloader <<<\n')

# collection notation: trailing slash or extension '.coll'
process_list = set([proc for proc in procs
                    if not (proc.endswith('/') or proc.endswith('.coll'))])
collections = [coll[:-1] + '.coll' for coll in procs if coll.endswith('/')]
collections.extend([coll for coll in procs if coll.endswith('.coll')])

# OpenLoops process API version
local_api_version = config['process_api_version']

repository_url = (config['remote_process_url'] +
                  '/%s/processes/' +  str(local_api_version))
collection_url = config['remote_process_url'] + '/%s/collections'
latest_api_version_url = (config['remote_process_url'] +
                          '/%s/processes/latest_version')
version_db_url = repository_url + '/versions.db'
channel_db_url = repository_url + '/channels.db'
libmappins_url = repository_url + '/librarymappings.db'
channel_db_file = os.path.join(config['process_lib_dir'], 'channels_%s.rinfo')


if sys.version_info[:2] < (2,5):
    def untar(archive, destinationpath):
        # There are still people around who use Python < 2.5
        #  which doesn't support the tarfile.extractall method.
        import subprocess
        if subprocess.call(['tar', '-xzf', archive, '--overwrite',
                            '--directory=' + destinationpath]) != 0:
            print('ERROR in untar: process archive extraction failed.')
            sys.exit(1)
else:
    def untar(archive, destinationpath):
        import tarfile
        lf = tarfile.open(archive, 'r')
        lf.extractall(destinationpath)
        lf.close()


def update_channel_db(repo):
    # get repository name of secret repository
    repo_name = OLToolbox.repo_name(repo)
    local_channel_file = channel_db_file % repo_name
    remote_channel_url = channel_db_url % repo
    if os.path.isfile(local_channel_file):
        fh = open(local_channel_file)
        local_hash = fh.readline().split()[0]
        fh.close()
    else:
        local_hash = None
    try:
        if remote_channel_url.startswith('/'):
            rfh = open(remote_channel_url, 'rb')
        else:
            rfh = urlopen(remote_channel_url)
    except (URLError, IOError) as e:
        print('Warning: Channel database update for repository ' + repo_name +
              ' failed (' + str(e) + '). Skip this repository.')
        return False
    hash_line = rfh.readline().decode()
    if local_hash != hash_line.split()[0]:
        local_hash = hashlib.md5()
        tmp_file = local_channel_file + '.~' + str(os.getpid())
        lfh = open(tmp_file, 'w')
        lfh.write(hash_line.strip() + '  ' +
                  time.strftime(OLToolbox.timeformat) + '\n')
        for line in rfh:
            lfh.write(line.decode())
            local_hash.update(line.strip())
        lfh.close()
        local_hash = local_hash.hexdigest()
        if local_hash == hash_line.split()[0]:
            os.rename(tmp_file, local_channel_file)
            print('updated channel database for repository', repo_name)
        else:
            print('ERROR: downloaded channel database inconsistent ' +
                  'for repository ' + repo_name)
            sys.exit(1)
    rfh.close()
    return True


def download(process, dbs, libmaps):
    # download process if one of the following conditions is met:
    # - the process API of the installed OpenLoops version
    #   differs from the one of the installed process;
    # - the version of the installed process (if any)
    #   differs from the version available on the server;
    # - download is forced by '--force' option;
    print('- process:', process, '...', end=' ')
    sys.stdout.flush()
    available = []
    for repo, cont in dbs.items():
        mapped_process = process
        if libmaps[repo] and process in libmaps[repo]:
            mapped_process = libmaps[repo][process]
            if mapped_process in downloaded.values():
                downloaded[process] = mapped_process
                print('mapped to', mapped_process, '(already downloaded)')
                return
            print('does not exist; downloading', mapped_process,
                  'instead ...', end=' ')
            sys.stdout.flush()
        if mapped_process in cont.content:
            # (process, date, hash, repo)
            available.append((mapped_process, cont.content[mapped_process][0],
                              cont.content[mapped_process][1], repo))
        elif mapped_process != process:
            print('\nOops! There seems to be something wrong with the ' +
                  'library mapping table. Please contact the authors.')
    # Order processes by date: latest last (to be used)
    available = sorted(
        available,
        key=lambda src: time.strptime(src[1], OLToolbox.timeformat))
    mapped_process = process
    if available:
        mapped_process = available[-1][0]
        downloaded[process] = mapped_process
    local_process_dir = os.path.join(config['process_src_dir'], mapped_process)
    version_installed = OLToolbox.import_dictionary(
        os.path.join(local_process_dir, 'version.info'), fatal = False)
    if version_installed is None:
        version_installed = {}
    try:
        api_installed = int(version_installed.get('process_api_version', None))
    except (AttributeError, TypeError):
        api_installed = None
    hash_installed = version_installed.get('hash', None)
    date_installed = version_installed.get('date', None)
    if not available:
        if args.ignore:
            print('IGNORED: not available (installed: %s)' % version_installed)
            return
        else:
            print('ERROR: not available (installed: %s)' % version_installed)
            sys.exit(1)
    # Only download, if the hash of the process to download differs
    # from the hash of the installed process, unless --force used.
    if not args.force:
        available = [src for src in available if src[2] != hash_installed]
    if not available:
        # the process is already available locally and is up-to-date.
        print('skipped: is up-to-date')
        return
    # Select the process which was uploaded last.
    available = available[-1]
    remote_archive = ((repository_url % available[3]) +
                      '/' + available[0] + '.tar.gz')
    local_archive = os.path.join(local_process_dir + '.tar.gz')
    # download the process
    print('download from repository: '+ available[3] +'...',end=' ')
    sys.stdout.flush()
    try:
        if remote_archive.startswith('/'):
            rf = open(remote_archive, 'rb')
        else:
            rf = urlopen(remote_archive)
    except (URLError, IOError):
        print('*** DOWNLOAD FAILED ***')
        if args.ignore:
            return
        else:
            sys.exit(1)
    lf = open(local_archive, 'wb')
    lf.write(rf.read())
    rf.close()
    lf.close()
    print('extract ...', end=' ')
    sys.stdout.flush()
    # remove target directory if it already exists
    try:
        shutil.rmtree(local_process_dir)
    except:
        pass
    # extract the process code from the archive
    untar(local_archive, config['process_src_dir'])
    # remove the local archive
    os.remove(local_archive)
    # store date and hash in the local process version.info
    process_version_file = os.path.join(local_process_dir, 'version.info')
    process_version = OLToolbox.import_dictionary(
        process_version_file, fatal=True)
    process_version['date'] = available[1]
    process_version['hash'] = available[2]
    OLToolbox.export_dictionary(process_version_file, process_version,
                                form = '%-25s %s')
    print('done')


process_dbs = {}
libname_maps = {}
max_latest_api_version = 0

if not os.path.isdir(config['process_lib_dir']):
    os.mkdir(config['process_lib_dir'])

first_repo = True
found_colls = set()
for repo in config['process_repositories']:
    repo_name = OLToolbox.repo_name(repo)
    # download the channel database for this repository
    if not update_channel_db(repo):
        continue
    # Latest OpenLoops process API version for which processes
    # are available in the repository
    latest_api_version = int(OLToolbox.import_dictionary(
        latest_api_version_url % repo)['process_api_version'])
    max_latest_api_version = max(max_latest_api_version, latest_api_version)
    # This fails if the repository does not contain processes
    # with the API of the installed OpenLoops version.
    process_dbs[repo] = OLToolbox.ProcessDB(db=(version_db_url % repo))
    libname_maps[repo] = OLToolbox.import_dictionary(libmappins_url % repo,
                                                     fatal=False)
    # scan all repositories for collections to download
    # Note that when the downloader is invoked by the build script,
    # the collections will already be resolved.
    for coll in collections:
        if coll == 'all.coll' or coll == repo_name + '.coll':
            process_coll = list(process_dbs[repo].content.keys())
        else:
            if first_repo:
                # check if the collection is available locally
                process_coll = OLToolbox.import_list(coll, fatal=False)
            else:
                process_coll = None
            if process_coll is None:
                # check if the collection is available in the repository
                process_coll = OLToolbox.import_list(
                    (collection_url % repo) + '/' + coll, fatal=False)
        if process_coll is not None:
            found_colls.add(coll)
        if process_coll:
            process_list.update(process_coll)
    first_repo = False

not_found_colls = [coll for coll in collections if coll not in found_colls]
if not_found_colls:
    if args.ignore:
        print('IGNORED: process collection(s) ', not_found_colls, ' not found.')
    else:
        print('ERROR: process collection(s) ', not_found_colls, ' not found.')
        sys.exit(1)

if local_api_version == max_latest_api_version:
    print('process API version: %d' % local_api_version)
elif local_api_version < max_latest_api_version:
    print('local process API version: %d (server: %d)' % (
          local_api_version, max_latest_api_version))
    print('******************************************')
    print('Please update your OpenLoops installation.')
    print('******************************************')

if not os.path.exists(config['process_src_dir']):
    os.mkdir(config['process_src_dir'])

downloaded = {}

for process in process_list:
    download(process, process_dbs, libname_maps)

OLToolbox.export_dictionary(
    os.path.join(config['process_src_dir'], 'downloaded.dat'), downloaded)

print('done\n')
