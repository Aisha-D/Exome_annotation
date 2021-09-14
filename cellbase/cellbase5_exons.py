from pandas.core.frame import DataFrame
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient
import pandas as pd
from pandas.core.frame import DataFrame
import urllib.request, json

from host import host_address


## 14/09/21
def flip_exonNumber(transcript_table):
    print("Negative strand! Flipping exons.")
    exonNumber = list(transcript_table.exonNumber)
    transcript_table2 = transcript_table.drop(columns=['exonNumber'])
    exonNumber.reverse()
    en = pd.DataFrame({"exonNumber" : exonNumber})
    transcript_table2 = pd.concat([transcript_table2, en], axis = 1)
    return(transcript_table2)


def query_cellbasedict(exon_dict):
    txs_dict = {}
    txs_dict['chr'] = exon_dict["chromosome"]
    txs_dict['exon_start'] = exon_dict["start"]
    txs_dict['exon_end'] = exon_dict["end"]
    txs_dict['gene_symbol'] = data['responses'][0]['results'][0]['name']
    txs_dict['transcript_id'] = data['responses'][0]['results'][0]['id'] #cannot take from exon as it has the _exonnumber attached to the transcript
    txs_dict['exonNumber'] = exon_dict["exonNumber"]
    #print(txs_dict)
    return txs_dict

# contains all transcripts
df = pd.DataFrame() 

g2t = pd.read_csv('g2t_210907_b38.tsv', sep='\t')

mane_transcripts = list(g2t.iloc[:,1])
len(mane_transcripts)
for transcript in range(10):
    print(mane_transcripts[transcript])
    url_get = "/webservices/rest/v5/hsapiens/feature/transcript/"
    txs = mane_transcripts[transcript]
    url_end = "/info?source=refseq"
    url_address = host_address + url_get + txs + url_end
    with urllib.request.urlopen(url_address) as url:
        data = json.loads(url.read().decode())
    
    #number of exons: 
    exons_num = len(data['responses'][0]['results'][0]['exons'])
    
    list_txs_dict = []
    for exon in range(exons_num):
        # print(exon)
        exon_dict = data['responses'][0]['results'][0]['exons'][exon]
        #print(exon_dict)
        if exon_dict["phase"] != -1:
            extracted_exons_dict = query_cellbasedict(exon_dict)
            list_txs_dict.append(extracted_exons_dict)
            # print(list_txs_dict)
        else:
            print("exon " + str(exon) + " is a non coding exon")
    
    # combine all exons into one table for that transcript
    transcript_table = pd.DataFrame(list_txs_dict)
    
    if data['responses'][0]['results'][0]['strand'] == "-":
        transcript_table2 = flip_exonNumber(transcript_table)
    
    df = df.append(transcript_table2) 

df.to_csv("exon_test.tsv", sep="\t", header=True, index=False)