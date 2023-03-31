#!/usr/bin/env python

from __future__ import print_function

# Collect all routines into a single file
# This is currently used for constructing the
# resulting combined libraries for siesta
# with LA

import argparse
import sys, os, re
import os.path as osp
import glob


def reduce_sort(lst):
    lst = list(set(lst))
    lst.sort()
    return lst
    

def list_read(files):
    lines = []
    for f in files:
        tmp = open(f).readlines()
        for line in tmp:
            line = line.strip()
            if len(line) == 0: continue
            line = line.replace('\n','')
            if not line.startswith('#'):
                lines.append(line.replace('\r',''))
    return lines


def strip_elements(lst, els):
    if els is None:
        return lst
    for el in els:
        try:
            lst.remove(el)
        except: pass
    return lst


# Luckily all files are in Fortran 77 format
# So it is easy to find subroutines used.
def calls_file(f):
    """ Reads file `f` and returns all subroutine calls from this file """

    # Only parse if it is a file
    if not osp.isfile(f): return []

    # Easy check for format
    is_f77 = f[-1] == 'f' or f[-1] == 'F'
    
    calls = []
    with open(f, 'r') as fh:
        for line in fh:
            
            # Remove all regular comments.
            line = line.lower().split('!')[0].strip()

            try:
                idx = line.index('call ')
            except:
                continue

            # We have a call in the line...
            call = line[idx+5:].split('(')[0]
            calls.append(call.strip())

    return calls

def calls_files(files):
    """ Reads all sources F, F90, f, f90 and returns a list of all calls
    made in all those files.
    """
    calls = []
    for f in files:
        calls.extend(calls_file(f))
    return calls


def calls_directory(d):
    """ Reads all sources F, F90, f, f90 and returns a list of all calls
    made in all those files.
    """
    # Only parse if it is a directory
    if not osp.isdir(d): return []

    calls = []
    for end in ['F', 'F90', 'f', 'f90']:
        files = glob.glob(osp.join(d,'*.'+end))
        calls.extend(calls_files(files))

    return calls


def subroutines_file(f):
    """ Reads file `f` and returns all subroutine definitions from this file """

    # Only parse if it is a file
    if not osp.isfile(f): return []
    
    subs = []
    with open(f, 'r') as fh:
        for line in fh:
            
            # Remove all regular comments.
            line = line.lower().split('!')[0].strip()

            try:
                idx = line.index('subroutine ')
                l = 10
            except:
                continue

            # We have a call in the line...
            try:
                sub = line[idx+l:].split('(')[0]
            except:
                continue
            subs.append(sub.strip())

    return subs


def find_files_calls(d, calls):
    """ Find files in `d` where the subroutines coinciding with those in `calls` may be found """

    # Only parse if it is a directory
    if not osp.isdir(d): return []

    files = []
    for end in ['F', 'F90', 'f', 'f90']:
        for f in glob.glob(osp.join(d,'*.'+end)):

            subs = subroutines_file(f)

            for sub in subs:
                if sub in calls:
                    files.append(f)
    return files


def find_files(d, infiles):
    """ Find files in `d` where the subroutines coinciding with those in `calls` may be found """

    # Only parse if it is a directory
    if not osp.isdir(d): return []

    files = []
    for infile in infiles:
        for f in glob.glob(osp.join(d,infile)):
            if osp.isfile(f):
                files.append(f)

    return files


def write_file(f, out, comments=False):
    with open(out, 'a') as oh:
        oh.write('! SOURCE-FILE = {}\n'.format(f))

        with open(f, 'r') as fh:
            for line in fh:
                if comments:
                    oh.write(line)
                elif line[0] in ' #':
                    oh.write(line)

def write_files(files, out, comments=False):
    for f in files:
        write_file(f, out, comments)
    
                    
if __name__ == '__main__':
    
    # First argument is the directory of the sources
    parser = argparse.ArgumentParser(description='Gather BLAS/LAPACK sources from external source directories.')

    # Add common options
    parser.add_argument('--directory','-d', dest="source_dir", type=str, nargs='+',
                        help='Specify directory of the sources where the entire source tree is found.')

    parser.add_argument('--file','-f', dest="source_file", type=str, nargs='+', default=[], 
                        help='Specify a specific file of sources.')

    parser.add_argument('--library-directory','-l', dest="library_dir", type=str, nargs='+', 
                        help='Specify directory of the sources where the entire source tree of the library is found.')

    parser.add_argument('--list-add-routine', dest='add_routine', type=str, nargs='+', default=[],
                        help='Read file with list of routines that are added as calls.')

    parser.add_argument('--list-remove-routine', dest='remove_routine', type=str, nargs='+', default=[],
                        help='Read file with list of routines that are removed as calls.')

    parser.add_argument('--list-add-file', dest='add_file', type=str, nargs='+', default=[],
                        help='Ensure these files are added to the source (found in [library_dir]).')

    parser.add_argument('--list-remove-file', dest='remove_file', type=str, nargs='+', default=[],
                        help='Ensure these files are removed from the source (found in [library_dir]).')


    parser.add_argument('--out','-o', dest='out_file', type=str,
                        help='Output file of resulting used files')

    parser.add_argument('--comments', action='store_true', default=False,
                        help='Strips all comments from the source files')

    args = parser.parse_args()

    # First remote out file
    try:
        os.remove(args.out_file)
    except: pass

    # Read contents of files on command-line
    add_routines = list_read(args.add_routine)
    rem_routines = list_read(args.remove_routine)
    add_files = list_read(args.add_file)
    rem_files = list_read(args.remove_file)


    # Retrieve all calls in the SOURCE directory
    calls = []
    for source_dir in args.source_dir:
        print('Parsing calls in: {}'.format(source_dir))
        calls.extend(calls_directory(source_dir))
    # Add custom files
    calls.extend(calls_files(args.source_file))

    # Add user-defined calls
    calls.extend(add_routines)
    calls = reduce_sort(calls)
    # Strip user-defined subroutines
    calls = strip_elements(calls, rem_routines)

    # Show how many different routines we actually found
    print('Found {} calls in the source directories...\n'.format(len(calls)))
    
    files = []
    for library_dir in args.library_dir:
        print('Parsing subroutines in: {}'.format(library_dir))
        files.extend(find_files_calls(library_dir, calls))

    # Now we have all higher level files directly
    # called from the source directories
    files = reduce_sort(files)
    len_source_files = len(files)

    # Add custom files
    tmp = []
    for library_dir in args.library_dir:
        files.extend(find_files(library_dir, add_files))
        tmp.extend(find_files(library_dir, rem_files))
    rem_files = reduce_sort(tmp)

    # Now keep finding routines untill all are found
    old_len = 0
    new_len = 1
    while old_len != new_len:
        old_len = len(files)

        # Be sure to have all sub-set files as well
        calls = calls_files(files)
        calls = reduce_sort(calls)
        # Strip user-defined calls
        calls = strip_elements(calls, rem_routines)
    
        print('Found {} calls in the library sources!'.format(len(calls)))
        for library_dir in args.library_dir:
            files.extend(find_files_calls(library_dir, calls))

        # Now we have all higher level files directly
        # called from the source directories
        # AND we have the library files directory calls
        files = reduce_sort(files)
        # Remove files
        for f in rem_files:
            try:
                files.remove(f)
            except:
                pass

        new_len = len(files)
    
    
    print('Require {} / {} files / + lib-files!'.format(len_source_files, len(files)))

    write_files(files, args.out_file, comments=args.comments)
    
