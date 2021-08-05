#!/usr/bin/python3
#
# 05/08/21


import pandas as pd
import argparse

def parse_args():
    """Reads on argument passed at the cmd line

    Returns:
        args: args held in an object
    """
    
    parser = argparse.ArgumentParser(description = 'bed files output from sql query')
    parser.add_argument(
        '-f', '--files',
        required=True,
        help='bed files'
        )

    args = parser.parse_args()

    return args

def read_bed(file):
    list_transcript_dict = []
    with open(file, 'r') as transcriptfile:
        for line in transcriptfile.readlines()[1:]:
            # The header for the file as follows:
            #['chrom', 'name', 'strand', 'txStart', 'txEnd', 'cdsStart', 
            # 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 
            # 'cdsStartStat', 'cdsEndStat', 'exonFrames']
            
            #make a dictionary with relevant headers for our final table
            transcript_dict = {}
            fields = line.strip().split("\t")
            transcript_dict['name'] = fields[1]
            transcript_dict['chrom'] = fields[0]
            transcript_dict['strand'] = fields[2]
            transcript_dict['exonCount'] = fields[7]
            transcript_dict['exonStarts'] = fields[8]
            transcript_dict['exonEnds'] = fields[9]
            transcript_dict['name2'] = fields[11]
            list_transcript_dict.append(transcript_dict)

    return list_transcript_dict

def make_bed(list_transcript_dict, df):

    # Loop through each transcript and filter the relevant
    # information. If the transcript is on the negative strand
    # reverse the order of the exon numbers.
    # At the end of each transcript, create a mini table with rows being
    # each exon of the transcript then append that mini table to the 
    # main dataframe with all the other transcripts
    for transcript_num in range(0, len(list_transcript_dict)):
        #print(transcript_num)
        tmp_dict = list_transcript_dict[transcript_num]
        lst = [[tmp_dict['name']]*int(tmp_dict['exonCount']), 
            [tmp_dict['chrom']]*int(tmp_dict['exonCount']),
            [tmp_dict['strand']]*int(tmp_dict['exonCount']),
            [tmp_dict['name2']]*int(tmp_dict['exonCount'])]
        transcript_table = pd.DataFrame(lst, index =['transcript', 'chr', 'strand', 'gene_symbol']).T
        exonStart =  tmp_dict['exonStarts'].split(",")[:-1]
        exonEnd = tmp_dict['exonEnds'].split(",")[:-1]
        exonNum = list(range(0, int(tmp_dict['exonCount']))) + [int(tmp_dict['exonCount'])]
        exonNum = exonNum[1:]
        if tmp_dict['strand'] == '-':
            #print("neg strand")
            exonNum.reverse()
            #print(exonNum)
        ES = pd.DataFrame({'exonStart': exonStart})
        ET = pd.DataFrame({'exonEnd': exonEnd})
        EN = pd.DataFrame({'exonNum': exonNum})
        transcript_table2 = pd.concat([transcript_table,ES,ET,EN], axis = 1)
        df = df.append(transcript_table2) 
        
    return df

def main():

    args = parse_args()
    input_filename = args.files
    list_transcript_dict = read_bed(input_filename)

    df = pd.DataFrame() 

    new_exome_df = make_bed(list_transcript_dict, df)

    # we are in the folder when running this we take the third part
    # remove the output.bed suffix of filename
    output_filename = input_filename.split("/")[2]
    output_filename = output_filename.split(".")[0]

    new_exome_df.to_csv(output_filename+'_reformatted.bed', sep="\t", header=True, index=False)

    print("Completed Reformatted file: " + output_filename+'_reformatted.bed')

if __name__ == "__main__":

    main()