#!/usr/bin/env python
#-*- coding: utf-8 -*-
#pylint: disable=
"""
File            : sequence_generation.py
Author          : Aurimas Repecka <aurimas.repecka AT gmail dot com>
Description     : Python script to generate sequences with minimal free energy.
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

APLHA = 0.5
TAU = 640

SEQUENCE_FILE = 'sequences/wt_p140.fasta'
OUTFILE_PATH = 'output/temp.fasta'
SVG_OUTFILE_PATH = 'output/'
STRUCTURE_PATH = 'secondary-structures/wt_p140.dat'
PARAMETERS_PATH = 'parameters/rna_andronescu2007.par'

NUCLEOTIDES_DICT = ['A', 'T', 'G', 'C']
NZONE_RANGE = [119,125]


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
    sorted_data = data[data[:,1].argsort()[::-1]]
    idxs = [i for i in range(len(sorted_data))]
    indexed_data = np.insert(sorted_data, 0, idxs, axis=1)

    trace0_data = indexed_data[indexed_data[:,3] == '0']
    trace1_data = indexed_data[indexed_data[:,3] == '1']
    trace2_data = indexed_data[indexed_data[:,3] == '2']
    trace3_data = indexed_data[indexed_data[:,3] == '3']
    trace4_data = indexed_data[indexed_data[:,3] == '4']
    trace5_data = indexed_data[indexed_data[:,3] == '5']
    trace6_data = indexed_data[indexed_data[:,3] == '6']

    trace0 = go.Scatter(
        x=trace0_data[:,0].astype(int),
        y=trace0_data[:,2].astype(float),
        mode='markers',
        name='0 mutations',
        text=trace0_data[:,1])
    trace1 = go.Scatter(
        x=trace1_data[:,0].astype(int),
        y=trace1_data[:,2].astype(float),
        mode='markers',
        name='1 mutations',
        text=trace1_data[:,1])
    trace2 = go.Scatter(
        x=trace2_data[:,0].astype(int),
        y=trace2_data[:,2].astype(float),
        mode='markers',
        name='2 mutations',
        text=trace2_data[:,1])
    trace3 = go.Scatter(
        x=trace3_data[:,0].astype(int),
        y=trace3_data[:,2].astype(float),
        mode='markers',
        name='3 mutations',
        text=trace3_data[:,1])
    trace4 = go.Scatter(
        x=trace4_data[:,0].astype(int),
        y=trace4_data[:,2].astype(float),
        mode='markers',
        name='4 mutations',
        text=trace4_data[:,1])
    trace5 = go.Scatter(
        x=trace5_data[:,0].astype(int),
        y=trace5_data[:,2].astype(float),
        mode='markers',
        name='5 mutations',
        text=trace5_data[:,1])
    trace6 = go.Scatter(
        x=trace6_data[:,0].astype(int),
        y=trace6_data[:,2].astype(float),
        mode='markers',
        name='6 mutations',
        text=trace6_data[:,1])
    data = [trace0, trace1, trace2, trace3, trace4, trace5, trace6]

    layout = go.Layout(  
        title = title,
        scene = go.Layout(
            xaxis = dict(title = xaxis),
            yaxis = dict(title = yaxis)
        )
    )

    fig = go.Figure(data = data, layout = layout)
    py.iplot(fig, filename=filename)



def generate_sequences_fold(sequence_path, structure_path, parameters_path, alpha, tau, nzone_range) :
    '''
    Generates different sequences and checks accuracy of folded structure
    param sequence_path: file path of analysed sequence
    param structure_path: file path of dot-bracket structure of analysed sequence
    param parameters_path: file path of parameters file for CoFold
    param alpha: alpha parameter of CoFold execution
    param tau: tau parameter of CoFold execution
    param nzone_range: range of nucleotides sequence that should be generated
    returns: array of Vector3 - (generated_fragment, accuracy percent, mutations count)
    '''

    results = []
    sequence_data = read_file(sequence_path).split("\n")
    base_structure = read_file(structure_path)
    seq_list = list(sequence_data[1])
    sequence_header = sequence_data[0] + '\n'
    nzone_fragment = seq_list[nzone_range[0]:nzone_range[1]]

    gen_fragments = itertools.product(NUCLEOTIDES_DICT, repeat=6)

    count = 0
    for gen_fragment in gen_fragments:
        seq_list[nzone_range[0]:nzone_range[1]] = gen_fragment
        new_seq = ''.join(seq_list)
        write_file(OUTFILE_PATH, sequence_header + new_seq)
        cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, OUTFILE_PATH)).read()

        sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]
        corr_perc = check_diff(base_structure, sec_struc)
        mutations = check_diff_mut(nzone_fragment, ''.join(gen_fragment))

        RNA.svg_rna_plot(new_seq, sec_struc, "%s%.5f_%s.svg" % (SVG_OUTFILE_PATH, corr_perc, str(count)))

        count = count + 1
        print("Iteration: %d" % count)
        results.append([''.join(gen_fragment), corr_perc, mutations])

    return np.asarray(results)



def main():
    
    results = generate_sequences_fold(SEQUENCE_FILE, STRUCTURE_PATH, PARAMETERS_PATH, 0.5, 640, NZONE_RANGE)
    #plot_scatter_chart("Generated sequences analysis (experiment target) -TGTAGC", results, "Iteration", "Accuracy", "140-sequence-analysis")

    #results_df = pd.DataFrame(results)
    #results_df.to_csv('output/svg/generated_seq_accuracies_exp.csv', header=False, index=False)

if __name__ == '__main__':
    main()