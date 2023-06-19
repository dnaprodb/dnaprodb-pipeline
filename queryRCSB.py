#!/usr/bin/env python
import argparse

import os
import json
import requests
from string import Template
from pymongo import MongoClient

from dnaprodb_utils import C

ROOT_DIR = C["ROOT_DIR"]

### Query RCSB for new entries
search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
query = open(os.path.join(ROOT_DIR, "external/structure_query.tmpl")).read()
query = Template(query)
query = json.loads(query.substitute(release_date='2019-10-01'))

req = requests.post(search_url, json=query)
RESULTS = req.json()

### Download structures
root_download_path = os.path.join(ROOT_DIR, "external/PDB/pdb_entries")
download_url = "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{}/{}.cif.gz"
downloaded = open(os.path.join(ROOT_DIR, "external/PDB/pdb_entries", "newest_downloaded_releases.txt"), "w")
for entry in RESULTS["result_set"]:
    entry_id  = entry["identifier"].lower()
    path = os.path.join(root_download_path, entry_id[0], "{}.cif.gz".format(entry_id))
    
    # download entry
    div = entry_id[1:3]
    try:
        r = requests.get(download_url.format(div, entry_id))
        with open(path, "wb") as FH:
            FH.write(r.content)
    except:
        continue
    downloaded.write(path+'\n')
downloaded.close()


# Add data to database
# client = MongoClient()
# db = client[DB_NAME]
# collection = db[COLLECTION_NAME]
