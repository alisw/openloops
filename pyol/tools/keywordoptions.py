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


import sys
import shlex

try:
    strtype = basestring
except NameError:
    strtype = str


class KeywordOptionsError(Exception):
    pass

class KeywordOption:
    pass


def is_string_list(ls):
    return isinstance(ls, (list, tuple, set)) and all(map(lambda el: isinstance(el, strtype), ls))


class KeywordOptions:

    def __init__(self, strict=True):
        self._options = []
        self._strict = strict

    def add(self, key, default=None, converter=None, group=None, required=False):
        """
        Add a keyword option.

        Arguments:
          key -- a keyword or if 'group' is specified a list of keywords.
          default -- the default value if the keyword is not given.
          converter -- a function to apply to the value before assigning it to the key.
          group -- the name of a dictionary to collect the keys in 'key'.
          required -- True: raise an error if the keyword is not given;
                      if 'group' is specified, True means that all keys are required;
                      also a subset of 'key' can be given.
        """
        if group is None:
            if not isinstance(key, strtype):
                raise ValueError("'key' must be a string, but is " + repr(key))
            if not isinstance(required, bool):
                raise ValueError("'required must be True or False, but is " + repr(required))
        elif isinstance(group, strtype):
            if default is None:
                default = {}
            if not is_string_list(key):
                raise ValueError("'key' must be a list of strings, but is " + repr(key))
            if not isinstance(required, bool) and not (is_string_list(required) and not (set(required) - set(key))):
                raise ValueError("'required must be True or False or a subset of key but is " + repr(required))
        else:
            raise ValueError("'group' must be None or a string, but is " + repr(group))


        option = KeywordOption()
        option.key = key
        option.default = default
        option.converter = converter
        option.group = group
        option.required = required
        self._options.append(option)

    def parse(self, args=None):
        """Parse arguments and attach keys/values to the KeywordOptions instance.

        'args' defaults to sys.argv. Also a string or list of strings can be given.
        If the arguments contain key=value, set self.key=value if key is not a group option.
        For group options set self.group[key]=value.
        'remaining_args' contains the list of non-keyword arguments
        and unrecognised arguments if strict=False.
        """
        if args is None:
            args = sys.argv

        elif isinstance(args, strtype):
            args = shlex.split(args)

        if not is_string_list(args):
            raise ValueError("Invalid arguments: " + str(args))

        splitargs = [arg.split("=", 1) for arg in args]
        args = dict([arg for arg in splitargs if len(arg) == 2])
        self.remaining_args = [arg[0] for arg in splitargs if len(arg) == 1]

        for opt in self._options:

            if isinstance(opt.key, strtype):
                if opt.key in args.keys():
                    val = args[opt.key]
                    del args[opt.key]
                    if opt.converter is not None:
                        val = opt.converter(val)
                else:
                    if opt.required:
                        raise KeywordOptionsError("Keyword option " + opt.key + " is required.")
                    else:
                        val = opt.default
                vars(self)[opt.key] = val
            else:
                # group option
                group = {}
                for key in opt.key:
                    if key in args.keys():
                        val = args[key]
                        del args[opt.key]
                        if opt.converter is not None:
                            val = opt.converter(val)
                        group[key] = val
                    else:
                        if (isinstance(opt.required, bool) and opt.required) or (
                              not isinstance(opt.required, bool) and key in opt.required):
                            raise KeywordOptionsError("Keyword option " + key + " is required.")
                if not group:
                    group = opt.default
                vars(self)[opt.group] = group

        self.unknown_options = args

        if self.unknown_options:
            if self._strict:
                raise KeywordOptionsError("Unrecognised arguments: " + ", ".join(unknown_options))
