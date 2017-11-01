#!/usr/bin/env python
#-*- coding: utf-8 -*-
#pylint: disable=
"""
File            : iterative_model.py
Author          : Aurimas Repecka <aurimas.repecka AT gmail dot com>
Description     : Python script to generate rna secondary structure using RNAfold library.
Packages        : https://www.tbi.univie.ac.at/RNA/#download 
"""

import os
import argparse
from datetime import datetime

from IPython import embed

import RNA

class OptionParser():
    def __init__(self):
        "User based option parser"
        self.parser = argparse.ArgumentParser(prog='iterative_model')
        self.parser.add_argument("--sequence", action="store", required=True,
            dest="sequence", default="", help="Destination of sequence file")
        self.parser.add_argument("--output", action="store", required=True,
            dest="output", default="", help="Directory where output will be stored")
        self.parser.add_argument("--window", action="store", type=int, required=True,
            dest="window", default=20, help="How much of sequence is shown in one iteration")
        self.parser.add_argument("--fixation", action="store", type=int, required=True,
            dest="fixation", default=10, help="Count of nucleotides that should be fixed after iteration")

def file_exists(file_path):
    """
    Checks if file exists
    :param file_path: path to file
    :raises: file not found error
    """

    if(not(os.path.isfile(file_path))):
        raise Exception("Sequence file in path %s was not found" % file_path)

def directory_exists(dir_path):
    """
    Checks if directory exists
    :param dir_path: path to directory
    :raises: directory not found error
    """

    if((not(os.path.exists(dir_path)))or(not(os.path.isdir(dir_path)))):
        raise Exception("Output directory in path %s was not found" % dir_path)

def read_file(file_path):
    """
    Reads given file into string
    :param file_path: path to file
    :raises: file not found error
    :returns: file content string
    """

    with open(file_path) as fin:
        fstr = fin.read()

    return fstr

def write_file(file_path, output):
    """
    Reads given file into string
    :param file_path: path to file
    :raises: file not found error
    """

    with open(file_path, 'w+') as fout:
        fout.write(output)

def init(ouput_dir, opts):
    """
    Initialization steps of program
    :param ouput_dir: path to output directory
    :param opts: program parameters
    :raises: directory not found error
    """

    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    working_dir = ouput_dir + "/" + now
    # Creates execution directory
    os.mkdir(working_dir)
    os.chdir(working_dir)

    # Writes parameters of execution to file
    write_file("parameters.txt", str(vars(opts)))

def write_iteration_results(seq, seq_count, rna_fold, final=False):
    """
    Writes iteration results: dot-bracket form, mfe and svg picture
    :param seq: string of sequence that was folded
    :param seq_count: integer that defines iteration count
    :param rna_fold: rna_fold objects
    """
    
    # Creates iteration directory
    os.mkdir(str(seq_count))

    # Creates dot-bracket and svg files
    output_files = "%d/%6.2f." % (seq_count, rna_fold[1])
    if(final):
        output_files = "%6.2f." % rna_fold[1]
    write_file(output_files + "dat", rna_fold[0])
    RNA.svg_rna_plot(seq, rna_fold[0], output_files + "svg")

def plot_results(range_, values):
    """
    Plots data using matplotlib
    :param range_: array of bins
    :param values: array of values
    """

    #TODO
    pass

#######################################################################################################

def main():
    "Main function"

    optpar  = OptionParser()
    opts = optpar.parser.parse_args()

    print("User input validation...")
    file_exists(opts.sequence)
    directory_exists(opts.output)

    print("Sequence reading...")
    seq = read_file(opts.sequence)

    print("Initialization...")
    init(opts.output, opts)

    print("Starting iteration process...")
    seq_len = len(seq)
    seq_start = 0
    seq_end = opts.window

    embed()

    '''
    while(seq_start < seq_len):
        seq_part = seq[0:seq_end]
        fc = RNA.fold_compound(seq_part)
        rna_fold = fc.mfe()

        write_iteration_results(seq_part, seq_start, rna_fold)
        print("-->%d:%d - %6.2f" % (seq_start, seq_end, rna_fold[1]))

        seq_start = seq_start + opts.fixation
        seq_end = min(seq_start + opts.window, seq_len-1)

    write_iteration_results(seq_part, seq_end, rna_fold, final=True)
    '''

if __name__ == '__main__':
    main()