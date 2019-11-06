//!******************************************************************************!
//! Copyright (C) 2014-2019 OpenLoops Collaboration. For authors see authors.txt !
//!                                                                              !
//! This file is part of OpenLoops.                                              !
//!                                                                              !
//! OpenLoops is free software: you can redistribute it and/or modify            !
//! it under the terms of the GNU General Public License as published by         !
//! the Free Software Foundation, either version 3 of the License, or            !
//! (at your option) any later version.                                          !
//!                                                                              !
//! OpenLoops is distributed in the hope that it will be useful,                 !
//! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
//! GNU General Public License for more details.                                 !
//!                                                                              !
//! You should have received a copy of the GNU General Public License            !
//! along with OpenLoops.  If not, see http://www.gnu.org/licenses/.             !
//!******************************************************************************!

#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <dlfcn.h>

#include <unistd.h>
#include <fcntl.h>

/* _______________ */
/* dirent wrappers */

// modes for dlopen (direct binding to Fortran doesn't work)
int ol_c_rtld_lazy = RTLD_LAZY;
int ol_c_rtld_now = RTLD_NOW;
int ol_c_rtld_global = RTLD_GLOBAL;
int ol_c_rtld_local = RTLD_LOCAL;

// global variable to hold a directory stream:
// only one directory can be open at a time.
DIR *ol_c_dirstream = NULL;

int ol_c_opendir(const char *dirname)
{
  if (ol_c_dirstream != NULL)
  {
    // printf("opendir failed: a directory is already open\n");
    return 127;
  };
  errno = 0;
  if ((ol_c_dirstream = opendir(dirname)) == NULL)
  {
    // printf("opendir: failed opening directory %s\n", dirname);
    return errno;
  };
  return 0;
}

int ol_c_readdir(char* entryname)
{
  struct dirent *direntry;
  entryname[0] = '\0';
  errno = 0;
  direntry = readdir(ol_c_dirstream);
  if (errno != 0)
  {
    // printf("readdir: reading directory content failed\n");
    return errno;
  };
  if (direntry != NULL)
  {
    strncpy(entryname, direntry->d_name, 256);
  };
  return 0;
}

void ol_c_closedir()
{
  closedir(ol_c_dirstream);
  ol_c_dirstream = NULL;
  return;
}

int ol_c_mkdir(const char *dirname)
{
  int err;
  err = mkdir(dirname, ACCESSPERMS);
  return err;
}

/* __________________ */
/* stdout redirection */

int ol_c_stdout_bak;

void ol_c_stdout_off()
{
  int devnull;
  fflush(stdout);
  ol_c_stdout_bak = dup(1);
  devnull = open("/dev/null", O_WRONLY);
  dup2(devnull, 1);
  close(devnull);
}

void ol_c_stdout_on()
{
  fflush(stdout);
  dup2(ol_c_stdout_bak, 1);
  close(ol_c_stdout_bak);
}
