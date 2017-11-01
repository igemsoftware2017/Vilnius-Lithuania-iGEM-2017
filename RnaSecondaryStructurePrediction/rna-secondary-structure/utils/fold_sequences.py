#!/usr/bin/env python
#-*- coding: utf-8 -*-
#pylint: disable=
"""
File            : fold_sequences.py
Author          : Aurimas Repecka <aurimas.repecka AT gmail dot com>
Description     : Python script to fold sequences with different lengths.
Packages        :   https://www.tbi.univie.ac.at/RNA/#download
                    http://www.e-rna.org/cofold/
"""

import os
import RNA
import argparse
import itertools
import numpy as np
import pandas as pd
from IPython import embed

import plotly.plotly as py
import plotly.graph_objs as go

SEQUENCE_FILE = '../../../sequences/'
TEMPFILE_PATH = '../../../temp/temp.fasta'
OUTFILES_PATH = '../_data/output/folded-mutations/'
STRUCTURE_PATH = '../../../secondary-structures/'
FOLDED_STRUCTURE_PATH = '../../../secondary-structures-folded/'
PARAMETERS_PATH = '../../../parameters/rna_andronescu2007.par'

STRUCTURE_LENGTHS = [132,140,160,200,555]
ALPHA = 0.5
TAU = 640


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


def check_diff(struc_1, struc_2):
    '''
    Check position-wise accuracy of two structures
    param struc_1: dot-bracket structure of sequence 1
    param struc_2: dot-bracket structure of sequence 2
    returns: similarity between two structures percentage
    '''
    corr = 0
    for idx, val in enumerate(struc_2):
        if val == struc_1[idx]:
            corr = corr + 1
    return corr/len(struc_1)

def check_diff_mut(seq_1, seq_2):
    '''
    Check position-wise accuracy of two sequences
    param struc_1: string of sequence 1
    param struc_2: string of sequence 2
    returns: mutations between two sequences
    '''

    mutations = 0
    for idx, val in enumerate(seq_1):
        if val != seq_2[idx]:
            mutations = mutations + 1
    return mutations

def parse_arguments():
    '''Parsing given arguments'''

    parser = argparse.ArgumentParser(description='Processing given mutations')
    parser.add_argument('--mutations', type=str, help='Mutations written in format: index1:mutation1_index2:mutation2...')
    args = parser.parse_args()
    
    mutations_str = args.mutations
    mutation_list = mutations_str.split('_')
    mutations = [mutation.split(':') for mutation in mutation_list]

    return mutations_str, mutations

def generate_sequences_fold(sequence_path, structure_path, parameters_path, alpha, tau, mutations_str,
    mutations, folded_structure_path) :
    '''
    Folds sequences with different mutations
    param sequence_path: file path of analysed sequence
    param structure_path: file path of dot-bracket structure of analysed sequence
    param parameters_path: file path of parameters file for CoFold
    param alpha: alpha parameter of CoFold execution
    param tau: tau parameter of CoFold execution
    param mutations_str: string of mutations in form of index:mutation_
    param mutations: list of mutations will in form of [[index, mutation]]
    param folded_structure_path: path to folded sequence dot-bracket structures
    returns: array of traces for visualization
    '''

    mutations_count = [0 for i in range(0, len(mutations))]
    os.chdir(OUTFILES_PATH)
    dirname = str(len(os.listdir())) + '_' + mutations_str
    os.mkdir(dirname)
    os.chdir(dirname)
    os.mkdir('sequences')
    os.mkdir('dot-bracket-structures')
    os.mkdir('folded-structures')

    for structure_length in STRUCTURE_LENGTHS:
        #Reading related sequence
        sequence_data = read_file(sequence_path + "wt_p" + str(structure_length) + ".fasta").split("\n")
        seq_list = list(sequence_data[1])
        sequence_header = sequence_data[0] + '\n'

        #Reading related structure
        base_structure = read_file(structure_path + "wt_p" + str(structure_length) + ".dat")
        base_structure_folded = read_file(folded_structure_path + "wt_p" + str(structure_length) + ".dat") 

        #Apply mutation
        mutated_sequence = seq_list
        for idx, mutation in enumerate(mutations):
            mutation_idx = int(mutation[0])
            mutation_seq = list(mutation[1])

            if mutation_idx > len(seq_list):
                continue
            elif mutation_idx + len(mutation_seq) > len(seq_list):
                mutation_seq  = mutation_seq[0:mutation_idx+len(mutation_seq)-len(seq_list)]

            base_seq = seq_list[mutation_idx:mutation_idx+len(mutation_seq)]
            mutations_count[idx] = check_diff_mut(base_seq, mutation_seq)
            mutated_sequence = seq_list
            mutated_sequence[mutation_idx:mutation_idx+len(mutation_seq)] = mutation_seq

        #Fold mutated sequence
        write_file(TEMPFILE_PATH, sequence_header + ''.join(mutated_sequence))
        cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, TEMPFILE_PATH)).read()
        sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]
        corr_perc = check_diff(base_structure, sec_struc)
        corr_perc_folded = check_diff(base_structure_folded, sec_struc)

        #Save results
        write_file('sequences/p%d.fasta' % (structure_length), '>' + mutations_str + '_p' + str(structure_length) + '\n' + ''.join(mutated_sequence))
        write_file('dot-bracket-structures/p%d.dat' % (structure_length), sec_struc)
        RNA.svg_rna_plot(''.join(mutated_sequence), sec_struc, "folded-structures/%d_%.5f_%.5f.svg" % (structure_length, corr_perc, corr_perc_folded))

    write_file('meta-data.txt', 'Index Mutation DiffCount \n' + str(list(zip(mutations,mutations_count))))


def main():
    mutations_str, mutations = parse_arguments()
    generate_sequences_fold(SEQUENCE_FILE, STRUCTURE_PATH, PARAMETERS_PATH, ALPHA, TAU, 
        mutations_str, mutations, FOLDED_STRUCTURE_PATH)
   
if __name__ == '__main__':
    main()