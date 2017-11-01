#!/usr/bin/env python
#-*- coding: utf-8 -*-
#pylint: disable=
"""
File            : cofold_analysis.py
Author          : Aurimas Repecka <aurimas.repecka AT gmail dot com>
Description     : Python script to analyse cofold parameters.
Packages        :   https://www.tbi.univie.ac.at/RNA/#download
                    http://www.e-rna.org/cofold/
"""

import os
import RNA
import numpy as np
from IPython import embed

import plotly.plotly as py
import plotly.graph_objs as go

APLHA_RANGE = 100
TAU_RANGE = 100

SEQUENCE_FILE = '../_data/output/sequence-generation/step-wise-accuracies/p132/wt_p132.fasta'
STRUCTURE_PATH = '../_data/output/sequence-generation/step-wise-accuracies/p132/wt_p132.dat'
PARAMETERS_PATH = '../_data/parameters/rna_andronescu2007.par'



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


def check_param_accuracy(alpha_range = APLHA_RANGE, tau_range = TAU_RANGE, 
                        sequence_file = SEQUENCE_FILE, structure_path = STRUCTURE_PATH, parameters_path = PARAMETERS_PATH):
    '''
    Loops between two parameters alpha and tau and checks cofold accuracy
    param alpha_range: number how many points between range 0 and 1 should be analysed as alpha
    param tau_range: number how many points between range 0 and 1000 should be analysed as tau
    param sequence_file: file path of analysed sequence
    param structure_path: path of correct dot-bracket structure of given sequence
    param parameters_path: path of parameters used by CoFold
    returns: arrays of alpha, tau and accuracy
    '''

    base_structure = read_file(structure_path)
    x = np.linspace(0,1,alpha_range)
    y = np.linspace(0,1000, tau_range)
    z, z_ = [], []

    for tau in y:
        for alpha in x:

            cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, sequence_file)).read()
            sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]
            corr_perc = check_diff(base_structure, sec_struc)
            print("%sx%s" % (str(alpha), str(tau)))
            z_.append(corr_perc)

        z.append(z_)
        z_ = []
    z = np.asarray(z)

    return x, y, z


def plot_surface_chart(title, x, y, z, xaxis, yaxis, zaxis, filename):
    '''
    Plots surface 3D chart in polt.ly servers
    param title: string title of chart name
    param x: data array of x axis
    param y: data array of y axis
    param z: data array of z axis
    param x: data array of x axis
    param xaxis: string name of x axis
    param yaxis: string name of x yaxis
    param zaxis: string name of x zaxis
    param filename: name of the file that will be saved in plot.ly
    '''

    data = go.Surface(x=x, y=y, z=z)
    layout = go.Layout(  
        title = title,
        scene = go.Scene(
            xaxis = dict(title = xaxis),
            yaxis = dict(title = yaxis),
            zaxis = dict(title = zaxis)
        )
    )

    fig = go.Figure(data = [data], layout = layout)
    py.iplot(fig, filename=filename)


def check_structure_accuracy(alpha, tau, sequence_file = SEQUENCE_FILE, structure_path = STRUCTURE_PATH, parameters_path = PARAMETERS_PATH):
    '''
    Checks sequence secondary structure accuracy and plots it in svg
    param alpha: alpha parameter of CoFold execution
    param tau: tau parameter of CoFold execution
    param sequence_file: file path of analysed sequence
    param structure_path: path of correct dot-bracket structure of given sequence
    param parameters_path: path of parameters used by CoFold
    '''

    base_structure = read_file(structure_path)
    cmd_out = os.popen('CoFold -d1 --noPS --distAlpha %.5f --distTau %.5f --paramFile=%s < %s' % (alpha, tau, parameters_path, sequence_file)).read()
    #cmd_out = os.popen('RNAfold -d1 --noPS < %s' % (sequence_file)).read()
    print(cmd_out)

    filename = cmd_out.split('\n')[0][1:]
    sequence = cmd_out.split('\n')[1]
    sec_struc = (cmd_out.split('\n')[2]).split(' ')[0]
    corr_perc = check_diff(base_structure, sec_struc)

    RNA.svg_rna_plot(sequence, sec_struc, "%s_%.5f.svg" % (filename, corr_perc))


def main():

    #x, y, z = check_param_accuracy(APLHA_RANGE, TAU_RANGE, SEQUENCE_FILE, STRUCTURE_PATH)
    #plot_surface_chart("CoFold Parameters analysis", x, y, z, 'Alpha', 'Tau', 'Accuracy(%)', 'CoFold Parameters analysis (Turner 1999)')

    check_structure_accuracy(0.5, 640, sequence_file = SEQUENCE_FILE, structure_path = STRUCTURE_PATH, parameters_path = PARAMETERS_PATH)

if __name__ == '__main__':
    main()