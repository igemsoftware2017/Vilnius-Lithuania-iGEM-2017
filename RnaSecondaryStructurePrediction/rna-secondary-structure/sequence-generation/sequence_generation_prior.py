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

SEQUENCE_FILE = '../_data/sequences/wt_p140.fasta'
SEQUENCE_FILE_C = '../_data/sequences/wt_p200.fasta'
SEQUENCE_FILE_C2 = '../_data/sequences/wt_p200.fasta'
TEMPFILE_PATH = '../_data/temp/temp.fasta'
SVG_OUTFILE_PATH = 'NAN'
STRUCTURE_PATH = '../_data/secondary-structures-folded/wt_p140.dat'
STRUCTURE_PATH_C = '../_data/secondary-structures-folded/wt_p200.dat'
STRUCTURE_PATH_C2 = '../_data/secondary-structures-folded/wt_p200.dat'
PARAMETERS_PATH = '../_data/parameters/rna_andronescu2007.par'

NUCLEOTIDES_DICT = ['A', 'T', 'G', 'C']
NUCLEOTIDES_SUBSTITUTION_DICT = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
NUMBERS_ARRAY = [0, 1, 2, 3, 4, 5, 6]

BASE_GROUP = [('T', 1), ('G', 2), ('T', 3), ('A', 4), ('G', 5), ('C', 6)]
PRIORITY_GROUP_1_O = [('G', 1), ('C', 2), ('G', 3), ('T', 4), ('T', 5), ('C', 6)]
PRIORITY_GROUP_2_O = [('C', 1), ('T', 2), ('C', 3), ('G', 4), ('G', 5), ('G', 6)]
PRIORITY_GROUP_3_O = [('A', 1), ('G', 2), ('A', 3), ('C', 4), ('C', 5), ('A', 6)]
PRIORITY_GROUP_4_O = [('T', 1), ('A', 2), ('T', 3), ('A', 4), ('A', 5), ('T', 6)]

# PRIORITIES CORRECT
PRIORITY_GROUP_1 = [('G', 1), ('A', 2), ('A', 3), ('C', 4), ('G', 5), ('C', 6)]
PRIORITY_GROUP_2 = [('C', 1), ('C', 2), ('C', 3), ('G', 4), ('A', 5), ('G', 6)]
PRIORITY_GROUP_3 = [('T', 1), ('G', 2), ('G', 3), ('T', 4), ('C', 5), ('T', 6)]
PRIORITY_GROUP_4 = [('A', 1), ('T', 2), ('T', 3), ('A', 4), ('T', 5), ('A', 6)]

CUGU_RANGE = [107,111]
NZONE_RANGE = [119,125]
NZONE_COMP_RANGE = [134,138]
NZONE_COMP_RANGE_2 = [80,86]

ALPHA = 0.5
TAU = 640
IS_COMPLEMENTARY = False
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
    traces = []
    sequence_data = read_file(sequence_path).split("\n")
    sequence_data_c2 = read_file(SEQUENCE_FILE_C2).split("\n")
    base_structure = read_file(structure_path)
    base_structure_c = read_file(STRUCTURE_PATH_C)
    base_structure_c2 = read_file(STRUCTURE_PATH_C2)
    seq_list = list(sequence_data[1])
    seq_list_c2 = list(sequence_data_c2[1])
    sequence_header = sequence_data[0] + '\n'
    nzone_fragment = seq_list[nzone_range[0]:nzone_range[1]]

    gen_praportions = itertools.product(NUMBERS_ARRAY, repeat=3)

    for praportion in gen_praportions:
        if sum(praportion) <= 1:
            data = []
            priority_1_count = praportion[0]
            priority_2_count = praportion[1]
            priority_3_count = praportion[2]
            folder_path = SVG_OUTFILE_PATH + "p" + str(praportion[0]) +  "-p" + str(praportion[1]) +  "-p" + str(praportion[2])
            if(not(IS_COMPARE)):
                os.mkdir(folder_path)
            gen_priority_1 = itertools.combinations(PRIORITY_GROUP_1, priority_1_count)

            for priority_1_group in gen_priority_1:
                priority_1_group = np.asarray(priority_1_group) if priority_1_group else np.empty(shape=(0, 2))
                priority_2_groups = [priority_2_element for priority_2_element in PRIORITY_GROUP_2 if str(priority_2_element[1]) not in priority_1_group[:,1]]
                gen_priority_2 = itertools.combinations(priority_2_groups, priority_2_count)

                for priority_2_group in gen_priority_2:
                    priority_2_group = np.asarray(priority_2_group) if priority_2_group else np.empty(shape=(0, 2))
                    priority_1_2_group = np.concatenate((priority_1_group[:,1], priority_2_group[:,1]))
                    priority_3_groups = [priority_3_element for priority_3_element in PRIORITY_GROUP_3 if str(priority_3_element[1]) not in priority_1_2_group]
                    gen_priority_3 = itertools.combinations(priority_3_groups, priority_3_count)

                    for priority_3_group in gen_priority_3:
                        priority_3_group = np.asarray(priority_3_group) if priority_3_group else np.empty(shape=(0, 2))
                        priority_1_2_3_group = np.concatenate((priority_1_group[:,1], priority_2_group[:,1], priority_3_group[:,1]))
                        priority_4_groups = [priority_4_element for priority_4_element in PRIORITY_GROUP_4 if str(priority_4_element[1]) not in priority_1_2_3_group]
                        gen_priority_4 = itertools.combinations(priority_4_groups, 6 - priority_3_count - priority_2_count - priority_1_count)

                        for priority_4_group in gen_priority_4:
                            priority_4_group = np.asarray(priority_4_group) if priority_4_group else np.empty(shape=(0, 2))

                            gen_fragment = np.concatenate((priority_1_group, priority_2_group, priority_3_group, priority_4_group))
                            sorted_fragment = gen_fragment[gen_fragment[:,1].argsort()[::1]][:,0]

                            seq_list[nzone_range[0]:nzone_range[1]] = sorted_fragment
                            seq_list_c2[nzone_range[0]:nzone_range[1]] = sorted_fragment

                            #if(is_complementary):
                            #    seq_list[NZONE_COMP_RANGE[0]:NZONE_COMP_RANGE[0] + 2] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0] + 4:nzone_range[1]]][::-1]
                            #    seq_list[NZONE_COMP_RANGE[0] + 2:NZONE_COMP_RANGE[1]] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0]:nzone_range[0] +2]][::-1]
                            #    seq_list_c2[NZONE_COMP_RANGE[0]:NZONE_COMP_RANGE[0] + 2] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0] + 4:nzone_range[1]]][::-1]
                            #    seq_list_c2[NZONE_COMP_RANGE[0] + 2:NZONE_COMP_RANGE[1]] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0]:nzone_range[0] +2]][::-1]
                            #    #seq_list[NZONE_COMP_RANGE_2[0]:NZONE_COMP_RANGE_2[1]] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0]:nzone_range[1]]][::-1]
                            #    #seq_list_c2[NZONE_COMP_RANGE_2[0]:NZONE_COMP_RANGE_2[1]] = [NUCLEOTIDES_SUBSTITUTION_DICT[elem] for elem in  seq_list[nzone_range[0]:nzone_range[1]]][::-1]

                            new_seq = ''.join(seq_list)
                            write_file(TEMPFILE_PATH, sequence_header + new_seq)
                            cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, TEMPFILE_PATH)).read()

                            sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]
                            corr_perc = check_diff(base_structure, sec_struc)
                            cugu_connected = not(sec_struc[CUGU_RANGE[0]:CUGU_RANGE[1]] == '....')

                            count = count + 1
                            print("Iteration: %d" % count)

                            if(not(IS_COMPARE)):
                                RNA.svg_rna_plot(new_seq, sec_struc, "%s/%.5f_%s_%s.svg" % (folder_path, corr_perc, ''.join(sorted_fragment), str(int(cugu_connected))))
                                data.append([''.join(sorted_fragment) + "_" + str(int(cugu_connected)), corr_perc])
                            else:
                                new_seq_c = new_seq[:-8]
                                write_file(TEMPFILE_PATH, sequence_header + new_seq_c)
                                cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, TEMPFILE_PATH)).read()

                                sec_struc_c = (cmd_out.split('\n')[2]).split(' ')[0]
                                corr_perc_c = check_diff(base_structure_c, sec_struc_c)
                                cugu_connected = not(sec_struc_c[CUGU_RANGE[0]:CUGU_RANGE[1]] == '....')

                                #new_seq_c2 = ''.join(seq_list_c2)
                                #write_file(TEMPFILE_PATH, sequence_header + new_seq_c2)
                                #cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, TEMPFILE_PATH)).read()

                                #sec_struc_c2 = (cmd_out.split('\n')[2]).split(' ')[0]
                                #corr_perc_c2 = check_diff(base_structure_c2, sec_struc_c2)
                                #cugu_connected = not(sec_struc_c2[CUGU_RANGE[0]:CUGU_RANGE[1]] == '....')
                                #data.append([corr_perc, ''.join(sorted_fragment) + "_" + str(int(cugu_connected)), corr_perc_c, corr_perc_c2])
                                data.append([corr_perc, ''.join(sorted_fragment) + "_" + str(int(cugu_connected)), corr_perc_c])

            data = np.asarray(data)
            sort_index = 0 if IS_COMPARE else 1
            sorted_data = data[data[:,sort_index].argsort()[::-1]]
            
            if(not(IS_COMPARE)):
                idxs = [i for i in range(count-len(sorted_data), count)]
                indexed_data = np.insert(sorted_data, 0, idxs, axis=1)
            else:
                indexed_data = sorted_data

            trace = go.Scatter(
                x=indexed_data[:,0].astype(float),
                y=indexed_data[:,2].astype(float),
                #z=indexed_data[:,3].astype(float),
                mode='markers',
                name="p" + str(praportion[0]) +  "-p" + str(praportion[1]) +  "-p" + str(praportion[2]),
                text=indexed_data[:,1])

            traces.append(trace)

    return traces


def main():
    
    traces = generate_sequences_fold(SEQUENCE_FILE, STRUCTURE_PATH, PARAMETERS_PATH, ALPHA, TAU, NZONE_RANGE, IS_COMPLEMENTARY)
    plot_scatter_chart("No 132", traces, "140 Accuracy", "132 Accuracy", "No 132")

    #results_df = pd.DataFrame(results)
    #results_df.to_csv('output/svg/generated_seq_accuracies_exp.csv', header=False, index=False)

if __name__ == '__main__':
    main()