#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 07:26:41 2019

@author: gmeier

mate generator
"""

from collections import defaultdict
import pysam


def read_pair_generator(bam, reference_name, start, end):
    """
    Find reads, store them in dict and return them once a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    # The reference parameter must match the header name of the sequence used for alignment.
    for read in bam.fetch(reference=reference_name, start=start, end=end):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
