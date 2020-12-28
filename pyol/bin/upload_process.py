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


#print('\n>>> OpenLoops Process Uploader <<<')

from __future__ import print_function
import sys

if sys.version_info[:2] < (2,7):
    print("This module requires Python 2.7 or later.")
    sys.exit(1)

import os
import shutil
import tarfile
import argparse
import hashlib
import time

sys.path.insert(0, os.path.abspath(os.path.join('pyol', 'config')))
sys.path.insert(0, os.path.abspath(os.path.join('pyol', 'tools')))

import OLBaseConfig
import OLToolbox

commandline_options = [arg.split('=',1) for arg in sys.argv[1:] if ('=' in arg and not arg.startswith('-'))]
config = OLBaseConfig.get_config(commandline_options)

if config['print_python_version']:
    print('upload_process.py uses Python', sys.version)

backup_old_processes = True

if not config['local_server_path']:
    print('option \'local_server_path\' must be specified in openloops.cfg')
    sys.exit(1)


parser = argparse.ArgumentParser(
    description="""Compress process code, upload it to a web server
                   and update the process version database.""")
parser.add_argument('processes', metavar='proc', nargs='*',
                    help="""process to upload, treated as a collection
                            if it ends with '/'""")
parser.add_argument('-d', '--delete', action='store_true',
                    help='delete processes from server')
parser.add_argument('-a', '--api', type=int, default=-1,
                    help='API version of process to delete (default: latest)')
parser.add_argument('-i', '--ignore', action='store_true',
                    help='ignore non-existing processes and collections')
parser.add_argument('-r', '--repository', metavar='repo', default=None,
                    help="""repository to which the processes will be added
                            (default: public)""")
parser.add_argument('-c', '--create', action='store_true',
                    help="""create a new repository, specified by '-r'.
                            If the repository name is followed by '*',
                            create a secret repository""")


def create_repository(repo, secret=False):
    if secret:
        import string
        import random
        seed_chars = string.ascii_lowercase + string.digits
        repo = (repo + '_' +
                ''.join(random.choice(seed_chars) for n in range(16)))
    repo_dir = os.path.join(config['local_server_path'], repo)
    process_dir = os.path.join(config['local_server_path'],
                               repo, 'processes')
    deprectated_dir = os.path.join(config['local_server_path'],
                                   repo, 'deprecated')
    collection_dir = os.path.join(config['local_server_path'],
                                  repo, 'collections')
    latest_api_version_file = os.path.join(process_dir, 'latest_version')
    if not os.path.isdir(repo_dir):
        os.mkdir(repo_dir)
    if not os.path.isdir(process_dir):
        os.mkdir(process_dir)
    if not os.path.isdir(deprectated_dir):
        os.mkdir(deprectated_dir)
    if not os.path.isdir(collection_dir):
        os.mkdir(collection_dir)
    OLToolbox.export_dictionary(latest_api_version_file,
                                {'process_api_version': 0})
    return repo


args = parser.parse_args(
    [arg for arg in sys.argv[1:] if (arg.startswith('-') or '=' not in arg)])

process_list = sum([proc.replace(',', ' ').split()
                    for proc in args.processes], [])
collections = [coll[:-1] for coll in process_list if coll.endswith('/')]
process_list = set([proc for proc in process_list if not proc.endswith('/')])

repository = args.repository
if not repository:
    repository = 'public'
if repository.endswith('*'):
    repository = repository[:-1]
    secret_repo = True
else:
    secret_repo = False

existing_repos = list(filter(
    lambda rep: os.path.isdir(os.path.join(config['local_server_path'], rep)),
    os.listdir(config['local_server_path'])))

if repository in existing_repos:
    repo_exists = True
    repo_key = repository
else:
    matching_secret = [rep for rep in existing_repos
                       if rep.startswith(repository + '_')]
    if len(matching_secret) == 0:
        repo_exists = False
    elif len(matching_secret) == 1:
        repo_exists = True
        repo_key = matching_secret[0]
    else:
        print('ambiguous repository name, matches')
        print(', '.join(matching_secret))
        sys.exit(1)

if args.create:
    if not args.repository:
        print('ERROR: repository to create must be explicitly specified.')
        sys.exit(1)
    if repo_exists:
        print('ERROR: repository already exists:', repo_key)
        sys.exit(1)
    else:
        repo_key = create_repository(repository, secret_repo)
        print('created repository', repo_key)
else:
    if not repo_exists:
        print('repository \'' + repository +
              '\' does not exist; use --create option to create it')
        sys.exit(1)


# repository structure, <local_server_path>/<repo>/
# directory for process archives/definitions which have been replaced
#   deprecated/
# directory for process collections
#   collections/
# processes, grouped in sub-directories by API version
#   processes/
#   latest process API version for which processes were uploaded
#     processes/latest_version
#   processes archives and definitions
#     processes/<api_version>/
#   process version database for uploaded processes
#     processes/<api_version>/version.db

# OpenLoops process API version
local_api_version = config['process_api_version']

repository_path = os.path.join(config['local_server_path'],
                               repo_key, 'processes', str(local_api_version))
deprecated_path = os.path.join(config['local_server_path'],
                               repo_key, 'deprecated')
collection_path = os.path.join(config['local_server_path'],
                               repo_key, 'collections')
latest_api_version_file = os.path.join(config['local_server_path'],
                                       repo_key, 'processes', 'latest_version')
# Latest OpenLoops process API version for which processes are available
# in the repository
latest_api_version = int(OLToolbox.import_dictionary(
    latest_api_version_file)['process_api_version'])
version_db_file = os.path.join(
    config['local_server_path'], repo_key, 'processes',
    str(local_api_version), 'versions.db')
channel_db_file = os.path.join(
    config['local_server_path'], repo_key, 'processes',
    str(local_api_version), 'channels.db')


def upload_process(process, db, ch_db, api):
    """Compress a process code folder and upload the archive
    to a repository on the web server."""
    # need: repository_path, deprecated_path, backup_old_processes
    print('- upload process:', process, '...', end=' ')
    sys.stdout.flush()
    old_date, old_hash, old_descr = db.get(process, (None, None, None))
    process_dir = os.path.join(config['process_src_dir'], process)
    process_version_file = os.path.join(process_dir, 'version.info')
    local_process_archive = process_dir + '.tar.gz'
    process_version = OLToolbox.import_dictionary(
        process_version_file, fatal=False)

    server_process_archive = os.path.join(repository_path, process + '.tar.gz')
    local_process_definition = os.path.join(config['process_def_dir'],
                                            process + '.m')
    server_process_definition = os.path.join(repository_path, process + '.m')
    server_backup_archive = os.path.join(
        deprecated_path, str(api) + '_' + process +
        '_' + str(old_date) + '_' + str(old_hash) + '.tar.gz')
    server_backup_definition = os.path.join(
        deprecated_path, str(api) + '_' + process +
        '_' + str(old_date) + '_' + str(old_hash) + '.m')

    if not process_version:
        if args.ignore:
            print('IGNORED: not available')
            return
        else:
            print('ERROR: not available')
            sys.exit(1)
    to_upload_api = int(process_version['process_api_version'])
    old_local_hash = process_version.get('hash', None)
    old_local_date = process_version.get('date', None)

    if to_upload_api != api:
        if args.ignore:
            print('IGNORED: process to upload does not match installed')
            print('         OpenLoops version (process: %d, OpenLoops: %d)' % (
                  to_upload_api, api))
            return
        else:
            print('ERROR: process to upload does not match installed')
            print('       OpenLoops version (process: %d, OpenLoops: %d)' % (
                  to_upload_api, api))
            sys.exit(1)

    if old_local_hash:
        # the local process was downloaded or uploaded before
        if old_local_hash == old_hash:
            print('skipped: is up-to-date')
            return
        elif old_date is not None and (
              time.strptime(old_local_date, OLToolbox.timeformat) <
              time.strptime(old_date, OLToolbox.timeformat)):
            print('skipped: process on server is newer')
            print('         (local: %s, server: %s)' % (
                  old_local_date, old_date))
            return

    if backup_old_processes and old_hash is not None:
        try:
            os.rename(server_process_archive, server_backup_archive)
        except OSError:
            print('[process backup failed]', end=' ')
            sys.stdout.flush()
        try:
            os.rename(server_process_definition, server_backup_definition)
        except OSError:
            print('[definition backup failed]', end=' ')
            sys.stdout.flush()

    # create process archive
    archive = tarfile.open(local_process_archive, 'w:gz')
    archive.add(process_dir, arcname=process)
    archive.close()
    # calculate archive hash and get upload time
    with open(local_process_archive, 'rb') as fh:
        archive_hash = hashlib.md5(fh.read()).hexdigest()
    upload_date = time.strftime(OLToolbox.timeformat)
    # store hash and upload time in local process directory
    process_version['date'] = upload_date
    process_version['hash'] = archive_hash
    OLToolbox.export_dictionary(process_version_file, process_version,
                                form = '%-25s %s')
    # get process description from process definition file
    with open(local_process_definition, 'r') as fh:
        description = OLToolbox.ProcessDB.no_description
        for line in fh:
            line = line.strip()
            if line.startswith('(*') and line.endswith('*)'):
                line = line[2:-2].strip()
                if line.startswith('description'):
                    description = line.split('=',1)[1].strip()
                    break
    # update process database
    db.update({process: (upload_date, archive_hash, description)})
    # update channel database
    info_options = OLToolbox.import_list(
        os.path.join(process_dir, 'info_' + process + '.txt'))
    info_options = [opt for opt in info_options if opt.startswith('options ')]
    if info_options:
        info_options = info_options[0].split()[1:]
    else:
        info_options = []
    info_files = OLToolbox.import_list(os.path.join(
        process_dir, "process_definition", "subprocesses.list"))
    info_files = [os.path.join(process_dir, "info_" + proc + ".txt")
                  for proc in info_files]
    info_files_extra = OLToolbox.import_list(os.path.join(
        process_dir, "process_definition", "subprocesses_extra.list"))
    info_files_extra = [os.path.join(process_dir, "info_" + proc + ".txt")
                  for proc in info_files_extra]
    channels = []
    for inf in info_files:
        channels.extend([line.split() + info_options
                         for line in OLToolbox.import_list(inf)])
    channels.sort(key=lambda el: el[1])
    channels_extra = []
    for inf in info_files_extra:
        channels_extra.extend([line.split() + info_options
                         for line in OLToolbox.import_list(inf)])
    channels_extra.sort(key=lambda el: el[1])
    ch_db.update({process: channels + channels_extra})
    # upload process archive and definition, delete temporary local archive
    shutil.copyfile(local_process_archive, server_process_archive)
    os.remove(local_process_archive)
    shutil.copyfile(local_process_definition, server_process_definition)
    print('done')


def delete_process(process, db, ch_db, api):
    """Delete a process from a repository on the server."""
    # need: repository_path, deprecated_path, backup_old_processes
    print('- delete process:', process, '...', end=' ')
    sys.stdout.flush()
    old_date, old_hash, old_descr = db.get(process, (None, None, None))
    if not old_date:
        print('skipped: does not exist')
        return
    server_process_archive = os.path.join(repository_path, process + '.tar.gz')
    server_process_definition = os.path.join(repository_path, process + '.m')
    server_backup_archive = os.path.join(
        deprecated_path, str(api) + '_' + process +
        '_' + str(old_date) + '_' + str(old_hash) + '.tar.gz')
    server_backup_definition = os.path.join(
        deprecated_path, str(api) + '_' + process +
        '_' + str(old_date) + '_' + str(old_hash) + '.m')
    if backup_old_processes:
        try:
            os.rename(server_process_archive, server_backup_archive)
        except OSError:
            print('[process backup failed]', end=' ')
            sys.stdout.flush()
        try:
            os.rename(server_process_definition, server_backup_definition)
        except OSError:
            print('[definition backup failed]', end=' ')
            sys.stdout.flush()
    else:
        try:
            os.remove(server_process_archive)
        except OSError:
            print('[deleting process failed]', end=' ')
            sys.stdout.flush()
            os.remove(server_process_definition)
        except OSError:
            print('[deleting definition failed]', end=' ')
            sys.stdout.flush()
    # update process and channel databases
    db.remove(process)
    ch_db.remove(process)
    print('done')


print('repository:', repository)

if local_api_version == latest_api_version:
    print('process API version: %d' % local_api_version)
    process_db = OLToolbox.ProcessDB(db=version_db_file)
    channel_db = OLToolbox.ChannelDB(db=channel_db_file)
else:
    if local_api_version > latest_api_version:
        print('new process API version: %d (server: %d)' % (
              local_api_version, latest_api_version))
        OLToolbox.export_dictionary(
            latest_api_version_file,
            {'process_api_version': local_api_version})
    else:
        print('WARNING: local process API is outdated')
        print('         (local: %d, server: %d)' % (
              local_api_version, latest_api_version))
    if not os.path.isdir(repository_path):
        os.mkdir(repository_path)
        process_db = OLToolbox.ProcessDB()
        channel_db = OLToolbox.ChannelDB()
        process_db.export_db(version_db_file)
        channel_db.export_db(channel_db_file)
    else:
        process_db = OLToolbox.ProcessDB(db=version_db_file)
        channel_db = OLToolbox.ChannelDB(db=channel_db_file)

for coll in collections:
    # Add content of collections to the process list.
    # The special collection 'all' selects all processes in the repository.
    if coll == 'all':
        process_coll = list(process_db.content.keys())
    else:
        process_coll = OLToolbox.import_list(
            os.path.join(collection_path, coll), fatal=False)
    if process_coll:
        process_list.update(process_coll)
    else:
        if args.ignore:
            print('IGNORED: process collection \'%s\' not found.' % coll)
        else:
            print('ERROR: process collection \'%s\' not found.' % coll)
            sys.exit(1)

for process in process_list:
    if args.delete:
        if args.api > 0:
            delete_api = args.api
        else:
            delete_api = latest_api_version
        delete_process(process, process_db, channel_db, delete_api)
    else:
        upload_process(process, process_db, channel_db, local_api_version)

if process_db.updated or channel_db.updated:
    print('update process database ...', end=' ')
    sys.stdout.flush()
    process_db.export_db(version_db_file)
    channel_db.export_db(channel_db_file)

print('done\n')
