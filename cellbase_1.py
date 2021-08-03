from pandas.core.frame import DataFrame
from pycellbase.cbclient import CellBaseClient
import pandas as pd
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description = 'Generate a bed file using cellbase annotation')
    parser.add_argument(
        '-g', '--genes',
        required=True,
        help='an array of genes (i.e. gene1, gene2 etc.)'
        )

    args = parser.parse_args()

    return args

def get_annotation(gc, genes, dataframe):

    for gene in genes:
        tf_response = gc.get_transcript(gene) 
        #tf_response_2 = tf_response[0]['result'] #select results section of tf infor
        #extract these following data from the gene transcript query
        transcript_chrom = []
        transcript_start = []
        transcript_end = []
        transcript_strand = []
        transcript_id = []
        transcript_name = []
        transcript_exon_number = []
        transcript_exon_id = []
        ## Loop through each transcript to get exon information
        for transcript in range(0,len(tf_response[0]['result'])):
            for exon in range(0,len(tf_response[0]['result'][transcript]['exons'])):
                transcript_id.append(tf_response[0]['result'][transcript]['id'])
                transcript_name.append(tf_response[0]['result'][transcript]['name'])
                transcript_exon_id.append(tf_response[0]['result'][transcript]['exons'][exon]['id'])
                transcript_exon_number.append(tf_response[0]['result'][transcript]['exons'][exon]['exonNumber'])
                transcript_chrom.append(tf_response[0]['result'][transcript]['exons'][exon]['chromosome'])
                transcript_start.append(tf_response[0]['result'][transcript]['exons'][exon]['start'])
                transcript_end.append(tf_response[0]['result'][transcript]['exons'][exon]['end'])
                transcript_strand.append(tf_response[0]['result'][transcript]['exons'][exon]['strand'])
        lst = [transcript_chrom, transcript_start, transcript_end, transcript_strand, 
                            transcript_id, transcript_name, transcript_exon_number, transcript_exon_id]
        gene_table = pd.DataFrame(lst, index =['chr', 'start', 'end', 'strand', 'transcript', 'transcript_name', 'exon_number', 'exon_id']).T
        dataframe = dataframe.append(gene_table) 

    return dataframe



def main():
    """
    Main function to generate bed file.
    """

    # To switch to version 5 as it has refseq annotations
    #custom_config = {'rest': {'hosts': ['bioinfo.hpc.cam.ac.uk:80/cellbase']}, 'version': 'v5', 'species': 'hsapiens'}
    #cc = ConfigClient(custom_config)
    #cbc = CellBaseClient(cc)
    #cbc.show_configuration()['version']

    # Default has v4
    cbc = CellBaseClient()
    gc = cbc.get_gene_client()

    dataframe = pd.DataFrame() 

    # read in genes
    #genes = ['BRCA1', 'APOE']
    args = parse_args()
    genes = list(filter(None, [x.strip() for x in args.genes.split(",")]))
    #print(list(args.genes))
    print("Genes: " + str(genes))

    annotated_genes = get_annotation(gc, genes, dataframe)

    filename = "_".join(genes) + '.bed'

    annotated_genes.to_csv(filename, sep="\t", header=True, index=False)

if __name__ == "__main__":

    main()