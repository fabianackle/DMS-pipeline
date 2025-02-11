# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 15:07:56 2018

@author: gmeier

DMS processing module

"""
import json
import Codon_truncation_bam_overlapp
import counter_fetch_single_position_paired_reads_in_progress
import create_count_file
import multiprocessing
from functools import partial
from contextlib import contextmanager


@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()


def DMS_processing(i, data_dict):
    print('analyzing files')
    Codon_truncation_bam_overlapp.codon_truncation(data_dict, i)
    counter_fetch_single_position_paired_reads_in_progress.count_mutants(data_dict, i)
#    create_count_file.make_HDF5(data_dict,i)


# running the dms data analysis using config data from jsonfile
def run_analysis(json_file_directory):
    with open(json_file_directory, 'r') as jsonfile:
        data_dict = json.load(jsonfile)

    with poolcontext(processes=15) as pool:
        pool.map(partial(DMS_processing, data_dict=data_dict), data_dict['data_files'])

    print('process finished')
