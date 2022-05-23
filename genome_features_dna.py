
import os, sys
import json
import re
import argparse
import json
from subprocess import Popen, PIPE, STDOUT
import time
import requests
import logging

def pretty_print_POST(req):
    """
    printed and may differ from the actual request.
    """
    print('{}\n{}\n{}\n\n{}'.format(
        '-----------START-----------',
        req.method + ' ' + req.url,
        '\n'.join('{}: {}'.format(k, v) for k, v in req.headers.items()),
        req.body,
    ))

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def genome_id_feature_gen(genome_ids, limit=250000):
    for gids in chunker(genome_ids, 4):
        selectors = ["eq(feature_type,rRNA)","eq(annotation,PATRIC)","in(genome_id,({}))".format(','.join(gids))]
        genomes = "and({})".format(','.join(selectors))   
        limit = "limit({})".format(limit)
        select = "sort(+genome_id,+sequence_id,+start)"
        base = "https://www.patricbrc.org/api/genome_feature/"
        query = "&".join([genomes, limit, select])
        headers = {"accept":"application/dna+fasta", "content-type": "application/rqlquery+x-www-form-urlencoded"}

        #Stream the request so that we don't have to load it all into memory
        r = requests.post(url=base, data=query, headers=headers, stream=True) 
        #r = requests.Request('POST', url=base, headers=headers, data=query)
        #prepared = r.prepare()
        #pretty_print_POST(prepared)
        #exit()
        if r.encoding is None:
            r.encoding = "utf-8"
        if not r.ok:
            logging.warning("Error in API request \n")
        for line in r.iter_lines(decode_unicode=True):
            yield line
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_id_files", type=str, nargs="*", default=["-"], help="Files with genome ids")
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()
    genome_ids=[]
    pargs = parser.parse_args()
    for f in pargs.genome_id_files:
        with open(f,'r') as fh:
            for line in fh:
                cur_ids=[i for i in re.split('; |, |,|\t|\n',line) if i]
                genome_ids.extend(cur_ids)
    for out_line in genome_id_feature_gen(genome_ids):
        print(out_line)

if __name__ == "__main__":
    start_time = time.time()
    main()
    sys.stderr.write(("--- %s seconds ---" % (time.time() - start_time))+"\n")