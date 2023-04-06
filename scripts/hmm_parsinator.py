#parsing Hmm tblout file
from collections import defaultdict
import pandas as pd
from Bio import SearchIO
import os,sys,argparse,subprocess
filename = sys.argv[1]
outname = sys.argv[2]

attribs = ['accession', 'bias', 'bitscore', 'description', 'cluster_num', 'domain_exp_num',  'domain_included_num', 'domain_obs_num', 'domain_reported_num', 'env_num', 'evalue', 'id', 'overlap_num', 'region_num']

hits = defaultdict(list)

with open(filename) as handle:
    for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
      #put query id in the dictionary
        hits['query_id'].append(queryresult.id)
      #print(queryresult.accession)
      #print(queryresult.description)
        for hit in queryresult.hits:
            for attrib in attribs:
                hits[attrib].append(getattr(hit, attrib))

pd.DataFrame.from_dict(hits).to_csv(outname, index=False)
