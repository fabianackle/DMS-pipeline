#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 21:06:22 2018

@author: gmeier
"""
import pysam
import os
import mate_generator
import mate_sequence_trimmer
# forward sequence trimmer function(reference sequence has to start with the first base in sequence ) include values for possible frame shift in sequence (start of frameshift and value for frame shift(1 or 2))
'''    
def sequence_trimmer_for(read,frameshift_position,frameshift_value=0):
    
    q=read.query_qualities
#store a shift value depending on the fact if the sequence start is bevore or after a frameshift in the sequence    
    if read.reference_start > frameshift_position:
        shift=frameshift_value-2
    else:
        shift=0
        
#no trimming needed if reference start 

    if read.reference_start in range(0+shift,4000,3):
        return read
       
    
    elif read.reference_start in range(2+shift,4000,3):
        read.query_sequence=read.query_sequence[1:]
        read.query_qualities=q[1:]
        read.reference_start=read.reference_start+1
        tuple_= read.cigartuples[0]
        length=tuple_[1]
        if length <10:
            read.cigarstring=str(length-1)+read.cigarstring[1:]
        elif length < 100:
            read.cigarstring=str(length-1)+read.cigarstring[2:]
        elif length < 1000:
            read.cigarstring=str(length-1)+read.cigarstring[3:]
        return read
             
        
    
    elif read.reference_start in range(1+shift,4000,3):
        read.query_sequence=read.query_sequence[2:]
        read.query_qualities=q[2:]
        read.reference_start=read.reference_start+2
        tuple_= read.cigartuples[0]
        length=tuple_[1]
        if length < 10:
            read.cigarstring=str(length-2)+read.cigarstring[1:]
        elif length < 100:
            read.cigarstring=str(length-2)+read.cigarstring[2:]
        elif length < 1000:
            read.cigarstring=str(length-2)+read.cigarstring[3:]
        return read



# function to trim query sequences on the right end
def sequence_trimmer_rev(read,frameshift_position,frameshift_value=0):
    q=read.query_qualities
    
    if read.reference_end > frameshift_position:
        shift=frameshift_value-2
    else:
        shift=0
        
        
    if read.reference_end in range(0+shift,4000,3):
        return read

    elif read.reference_end in range(2+shift,4000,3):
        read.query_sequence=read.query_sequence[:-2]
        read.query_qualities=q[:-2]
#        read.reference_end=read.reference_end+1
        tuple_= read.cigartuples[-1]
        length=tuple_[1]
        if length < 10:
            read.cigarstring=read.cigarstring[0:-2]+str(length-2)+read.cigarstring[-1]
        elif length < 100:
            read.cigarstring=read.cigarstring[0:-3]+str(length-2)+read.cigarstring[-1]
        elif length < 1000:
            read.cigarstring=read.cigarstring[0:-4]+str(length-2)+read.cigarstring[-1]
        return read
    
    
    elif read.reference_end in range(1+shift,4000,3):
        read.query_sequence=read.query_sequence[:-1]
        read.query_qualities=q[:-1]
#        read.reference_end=read.reference_end+2
        tuple_= read.cigartuples[-1]
        length=tuple_[1]
        if length < 10:
            read.cigarstring=read.cigarstring[0:-2]+str(length-1)+read.cigarstring[-1]
        elif length < 100:
            read.cigarstring=read.cigarstring[0:-3]+str(length-1)+read.cigarstring[-1]
        elif length < 1000:
            read.cigarstring=read.cigarstring[0:-4]+str(length-1)+read.cigarstring[-1]
        return read
'''


def codon_truncation(data_dict, curent_file):

    input_file_path = data_dict['inputdir'] + '/' + curent_file
    output_directory = data_dict['outputdir'] + '/' + curent_file[:-4] + '_codontruncated.bam'
    freameshift_position = data_dict['frameshift_position']
    frameshift_offset = data_dict['frameshift_offset']
# indexing inputfile
    os.system('samtools index ' + input_file_path.replace(' ', '\ ') + ' >' + input_file_path.replace(' ', '\ ') + '.bai')
    # open file for reading and new file for writing

    samfile = pysam.AlignmentFile(input_file_path, 'rb')
    samfile_trimmed = pysam.AlignmentFile(output_directory, 'wb', template=samfile)

    # remove all hardclipped reads and all unmapped reads

    unmapped_read_count = 0
    softclip_count = 0
    hardclip_count = 0
    short_overlap_or_non = 0
    for read1, read2 in mate_generator.read_pair_generator(samfile, start=None, end=None):

        if read1.flag == 4 or read2.flag == 4 or read1.cigarstring == None or read2.cigarstring == None:
            unmapped_read_count = unmapped_read_count + 1
#        elif 'S' in read1.cigarstring or 'S' in read2.cigarstring:
#            softclip_count=softclip_count+1
        elif 'H' in read1.cigarstring or 'H' in read2.cigarstring:
            hardclip_count = hardclip_count + 1
        else:
            #            val=False
         #           if read1.query_name=='A00460:222:HF5YKDRXX:1:2242:5701:25113':
            #            if 'S' in read1.cigarstring and 'D' in read1.cigarstring:
            #            if read1.reference_start!= read2.reference_start:
            #                val=True
            #                print('\n')
            #                print(read1)
            #                print(read1.cigartuples)
            #                print(read1.cigarstring)
            #                print(read2)
            #                print(read2.cigartuples)
            #                print(read2.cigarstring)
            trimmed_read1, trimmed_read2 = mate_sequence_trimmer.trim(read1, read2, freameshift_position, frameshift_offset)


#            if trimmed_read1 ==None:
#                print('damn')
#            forwardtrimmed_read=sequence_trimmer_for(read,freameshift_position,frameshift_offset)
#            for_rev_trimmed_read=sequence_trimmer_rev(forwardtrimmed_read,freameshift_position,frameshift_offset)
            if trimmed_read1 == None or trimmed_read2 == None:
                short_overlap_or_non = short_overlap_or_non + 1

            elif trimmed_read1 != None and trimmed_read2 != None:
                #                if read1.query_name=='A00460:222:HF5YKDRXX:1:2220:30110:8484':
                #                if 'I' in read1.cigarstring:
                #                if val==True:
                #                if trimmed_read1.reference_length!= trimmed_read2.reference_length:
                #                if trimmed_read1.query_length!=trimmed_read2.query_length:

                #                    print(trimmed_read1.query_length)
                #                    print(trimmed_read2.query_length)
                #                    print(trimmed_read1)
                #                    print(trimmed_read1.cigartuples)
                #                    print(trimmed_read1.cigarstring)
                #                    print(trimmed_read2)
                #                    print(trimmed_read2.cigartuples)
                #                    print(trimmed_read2.cigarstring)
                #                    print('\n')
                samfile_trimmed.write(trimmed_read1)
                samfile_trimmed.write(trimmed_read2)
    print('reads processed')

    samfile.close()
    samfile_trimmed.close()

    # indexing output file
    os.system('samtools index ' + output_directory.replace(' ', '\ ') + ' >' + output_directory.replace(' ', '\ ') + '.bai')

    print(str(unmapped_read_count) + 'unmapped reads, ' + str(hardclip_count) + ' hardclipped reads and ' + str(short_overlap_or_non) + ' short overlap reads were removed')
