#!/usr/bin/env python
#-*- coding: utf-8 -*-
#pylint: disable=
"""
File            : sequence_generation_prior.py
Author          : Aurimas Repecka <aurimas.repecka AT gmail dot com>
Description     : Python script to generate sequences with minimal free energy (priority groups).
Packages        :   https://www.tbi.univie.ac.at/RNA/#download
                    http://www.e-rna.org/cofold/
"""

import os
import RNA
import itertools
import numpy as np
import pandas as pd
from IPython import embed

import plotly.plotly as py
import plotly.graph_objs as go

SEQUENCE_FILE = '../_data/sequences/'
TEMPFILE_PATH = '../_data/temp/temp.fasta'
SVG_OUTFILE_PATH = '../_data/output/sequence-generation/1-length-deletion/complementary/'
STRUCTURE_PATH = '../_data/secondary-structures/'
PARAMETERS_PATH = '../_data/parameters/rna_andronescu2007.par'

NUCLEOTIDES_DICT = ['A', 'T', 'G', 'C']
NUCLEOTIDES_SUBSTITUTION_DICT = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
NUMBERS_ARRAY = [0, 1, 2, 3, 4, 5, 6]

BASE_GROUP = [('T', 1), ('G', 2), ('T', 3), ('A', 4), ('G', 5), ('C', 6)]
PRIORITY_GROUP_1 = [('G', 1), ('A', 2), ('A', 3), ('C', 4), ('G', 5), ('C', 6)]
PRIORITY_GROUP_2 = [('C', 1), ('C', 2), ('C', 3), ('G', 4), ('A', 5), ('G', 6)]
PRIORITY_GROUP_3 = [('T', 1), ('G', 2), ('G', 3), ('T', 4), ('C', 5), ('T', 6)]
PRIORITY_GROUP_4 = [('A', 1), ('T', 2), ('T', 3), ('A', 4), ('T', 5), ('A', 6)]


CUGU_RANGE = [107,111]
NZONE_RANGE = [119,125]
NZONE_COMP_RANGE = [134,138]
STRUCTURE_LENGTHS = [132,140,200,555]

ALPHA = 0.5
TAU = 640
IS_COMPLEMENTARY = True
IS_COMPARE = True


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


def plot_scatter_chart(title, data, xaxis, yaxis, filename):
    '''
    Plots scatter chart in polt.ly servers
    param title: string title of chart name
    param x: data array of x axis
    param y: data array of y axis
    param xaxis: string name of x axis
    param yaxis: string name of y yaxis
    param filename: name of the file that will be saved in plot.ly
    '''

    layout = go.Layout(  
        title = title,
        scene = go.Layout(
            xaxis = dict(title = xaxis),
            yaxis = dict(title = yaxis)
        )
    )

    fig = go.Figure(data = data, layout = layout)
    py.iplot(fig, filename=filename)


def read_sequences_list(sequence_path):
    '''
    Reads multiple sequences in the list
    param sequence_path: path of sequences locations
    returns list of sequences
    '''

    res = []
    for structure_length in STRUCTURE_LENGTHS:
        sequence_data = read_file(sequence_path + "wt_p" + str(structure_length) + ".fasta").split("\n")
        seq_list = list(sequence_data[1])
        sequence_header = sequence_data[0] + '\n'
        res.append((seq_list, sequence_header))
    return res




def generate_sequences_fold(sequence_path, structure_path, parameters_path, alpha, tau, nzone_range, is_complementary) :
    '''
    Generates different sequences and checks accuracy of folded structure
    param sequence_path: file path of analysed sequence
    param structure_path: file path of dot-bracket structure of analysed sequence
    param parameters_path: file path of parameters file for CoFold
    param alpha: alpha parameter of CoFold execution
    param tau: tau parameter of CoFold execution
    param nzone_range: range of nucleotides sequence that should be generated
    param is_complementary: parameter for complementary sequence adjustments
    param mutations: list of mutations will be applied in nzone range
    returns: array of traces for visualization
    '''

    traces = []
    sequence_data = read_sequences_list(sequence_path)
    base_structures = [read_file(structure_path + "wt_p" + str(structure_length_) + ".dat") for structure_length_ in STRUCTURE_LENGTHS]

    for i in range(NZONE_RANGE[0],NZONE_RANGE[1]):
        structure_res = []
        folder_path = SVG_OUTFILE_PATH + str(i)
        os.mkdir(folder_path)

        for structures_ in zip(STRUCTURE_LENGTHS, base_structures, sequence_data):
            structure_length = structures_[0]
            base_structure_ = structures_[1]
            seq_list, sequence_header = structures_[2]
            seq_list_ = list(seq_list)
            print("%s_%d" % (str(i), structure_length))

            seq_list_.pop(i)
            if(is_complementary):
                seq_list_.pop(i-34)

            new_seq = ''.join(seq_list_)
            write_file(TEMPFILE_PATH, sequence_header + new_seq)
            cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, TEMPFILE_PATH)).read()
            sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]

            corr_perc = check_diff(base_structure_, sec_struc)

            RNA.svg_rna_plot(new_seq, sec_struc, "%s/%d_%.5f_%s.svg" % (folder_path, structure_length, corr_perc, str(i)))
            structure_res.append([structure_length, corr_perc])

        data = np.asarray(structure_res)
        sorted_data = data[data[:,0].argsort()[::1]]
        trace = go.Scatter(
            x=sorted_data[:,0].astype(int),
            y=sorted_data[:,1].astype(float),
            mode='markers',
            name=str(i))
        traces.append(trace)

    return traces


def main():
    traces = generate_sequences_fold(SEQUENCE_FILE, STRUCTURE_PATH, PARAMETERS_PATH, ALPHA, TAU, NZONE_RANGE, IS_COMPLEMENTARY)
    plot_scatter_chart("1 length complementary deletion", traces, "Length of structure", "Accuracy", "1 length complementary deletions ")

if __name__ == '__main__':
    main()