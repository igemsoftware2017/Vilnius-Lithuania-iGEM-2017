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
SVG_OUTFILE_PATH = 'NAN'
STRUCTURE_PATH = '../_data/secondary-structures/'
PARAMETERS_PATH = '../_data/parameters/rna_andronescu2007.par'

NUCLEOTIDES_DICT = ['A', 'U', 'C', 'G']
NUCLEOTIDES_SUBSTITUTION_DICT = {'A' : 'U', 'C' : 'G', 'G' : 'C', 'U' : 'A'}
NUMBERS_ARRAY = [0, 1, 2, 3, 4, 5, 6]

NZONE_RANGE = [141,147]
NZONE_COMP_RANGE = [113,115]
STRUCTURE_LENGTHS = [160,200,555]

ALPHA = 0.5
TAU = 640
IS_COMPLEMENTARY = False


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
            yaxis = dict(title = yaxis),
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
    returns: array of traces for visualization
    '''

    count = 0
    res = []
    traces = []
    sequence_data = read_sequences_list(sequence_path)
    base_structures = [read_file(structure_path + "wt_p" + str(structure_length_) + ".dat") for structure_length_ in STRUCTURE_LENGTHS]
    structures__ = list(zip(STRUCTURE_LENGTHS, base_structures, sequence_data))
    nzone_fragment = sequence_data[0][0][nzone_range[0]:nzone_range[1]]

    gen_fragments = itertools.product(NUCLEOTIDES_DICT, repeat=6)

    for fragment in gen_fragments:
        mutations = check_diff_mut(nzone_fragment, ''.join(fragment))
        structure_res = [''.join(fragment), mutations]

        for structures_ in structures__:
            structure_length = structures_[0]
            base_structure_ = structures_[1]
            seq_list, sequence_header = structures_[2]

            seq_list[nzone_range[0]:nzone_range[1]] = fragment

            if(is_complementary):
                seq_list[NZONE_COMP_RANGE[0]:NZONE_COMP_RANGE[1]] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0] + 4:nzone_range[1]]][::-1]

            new_seq = ''.join(seq_list)
            write_file(TEMPFILE_PATH, sequence_header + new_seq)
            cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, TEMPFILE_PATH)).read()
            sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]
            corr_perc = check_diff(base_structure_, sec_struc)

            structure_res.append(corr_perc)

        res.append(structure_res)
        print("Iteration: %d Results: %s" % (count, str(structure_res)))
        count = count + 1

    data = np.asarray(res)
    trace0_data = data[data[:,1] == '0']
    trace1_data = data[data[:,1] == '1']
    trace2_data = data[data[:,1] == '2']
    trace3_data = data[data[:,1] == '3']
    trace4_data = data[data[:,1] == '4']
    trace5_data = data[data[:,1] == '5']
    trace6_data = data[data[:,1] == '6']

    trace0 = go.Scatter3d(
        x=trace0_data[:,2].astype(float),
        y=trace0_data[:,3].astype(float),
        z=trace0_data[:,4].astype(float),
        mode='markers',
        name='0 mutations',
        text=trace0_data[:,0])
    trace1 = go.Scatter3d(
        x=trace1_data[:,2].astype(float),
        y=trace1_data[:,3].astype(float),
        z=trace1_data[:,4].astype(float),
        mode='markers',
        name='1 mutations',
        text=trace1_data[:,0])
    trace2 = go.Scatter3d(
        x=trace2_data[:,2].astype(float),
        y=trace2_data[:,3].astype(float),
        z=trace2_data[:,4].astype(float),
        mode='markers',
        name='2 mutations',
        text=trace2_data[:,0])
    trace3 = go.Scatter3d(
        x=trace3_data[:,2].astype(float),
        y=trace3_data[:,3].astype(float),
        z=trace3_data[:,4].astype(float),
        mode='markers',
        name='3 mutations',
        text=trace3_data[:,0])
    trace4 = go.Scatter3d(
        x=trace4_data[:,2].astype(float),
        y=trace4_data[:,3].astype(float),
        z=trace4_data[:,4].astype(float),
        mode='markers',
        name='4 mutations',
        text=trace4_data[:,0])
    trace5 = go.Scatter3d(
        x=trace5_data[:,2].astype(float),
        y=trace5_data[:,3].astype(float),
        z=trace5_data[:,4].astype(float),
        mode='markers',
        name='5 mutations',
        text=trace5_data[:,0])
    trace6 = go.Scatter3d(
        x=trace6_data[:,2].astype(float),
        y=trace6_data[:,3].astype(float),
        z=trace6_data[:,4].astype(float),
        mode='markers',
        name='6 mutations',
        text=trace6_data[:,0])
    traces = [trace0, trace1, trace2, trace3, trace4, trace5, trace6]

    return traces


def main():
    
    traces = generate_sequences_fold(SEQUENCE_FILE, STRUCTURE_PATH, PARAMETERS_PATH, ALPHA, TAU, NZONE_RANGE, IS_COMPLEMENTARY)
    plot_scatter_chart("-35 region - one-sided -new", traces, "200 Accuracy", "555 Accuracy", "-35 region - one-sided -new")

    #results_df = pd.DataFrame(results)
    #results_df.to_csv('output/svg/generated_seq_accuracies_exp.csv', header=False, index=False)

if __name__ == '__main__':
    main()