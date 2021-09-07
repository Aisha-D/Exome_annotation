from pandas.core.frame import DataFrame
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient
import pandas as pd
import argparse
import urllib.request, json 

from host import host_address

## Query cell database it diffucult therefore all querying will be
## done via url which will take time but ensure its the correct version 5


# Read in the HGNC list that currently exists 
HGNC_df = pd.read_csv('HGNC_210902.tsv', sep='\t')
# this list will contain all the HGNC_ids and relevant MANE RefSeq transcripts
all_genes = []

# len(HGNC_df)

for row in range(10):
    HGNC_id = HGNC_df.iloc[row, 0] #HGNC ID is in 0th column
    ensemble_id = HGNC_df.iloc[row, 11] #Ensembl gene ID is in 12th column
    # make dict to keep info on what refseq is for ensemble gene in cellbase
    gene_dict = {}
    gene_dict['HGNC_ID'] = HGNC_id
    print(HGNC_id)
    ## get url version
    url_get = "feature/gene/"
    url_end = "/info"
    url_address = host_address + url_get + ensemble_id + url_end
    with urllib.request.urlopen(url_address) as url:
        data = json.loads(url.read().decode())
    # how many transcripts are there?
    txs_num = len(data['responses'][0]['results'][0]['transcripts'])
    # for loop the ensemble transcripts to find what the refseq mane 
    # transcripts are for the ensemble gene
    for transcript in range(txs_num):
        dicts = list(data['responses'][0]['results'][0]['transcripts'][transcript]['xrefs'])
        mane_dict = [item for item in dicts if item["dbName"] == "mane_select_refseq"]
        if not mane_dict:
            # some ensemble transcripts do not have refseq mane transcript
            # for these, we skip to the next transcript
            continue
        else:
            # print("Transcript number " + str(transcript) + " DOES have Mane RefSeq Annotaion")
            gene_dict['MANE_RefSeqID'] = mane_dict[0]['id']
    # append the gene_dict to the list of all HGNC ids
    all_genes.append(gene_dict) 



df = pd.DataFrame.from_dict(all_genes)

    
