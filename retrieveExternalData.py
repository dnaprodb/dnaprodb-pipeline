#!/usr/bin/env python
import os
import json
import argparse
import requests
import pickle
#import re
#import xmltodict
from string import Template
from dnaprodb_utils import C
#from Bio.PDB.MMCIF2Dict import MMCIF2Dict
#from Bio import SwissProt

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("pdb_list_file")
arg_parser.add_argument("--uniprot_file", help="a JSON file of previously downloaded UniProt entries")
arg_parser.add_argument("--no_refresh_uniprot", dest="refresh_uniprot", action='store_false', help="Do not update UniProt file with newly downloaded data")
arg_parser.add_argument("--no_download_clusters", dest="download_clusters", action='store_false', help="Do not update cluster mappings with newly downloaded data")
ARGS = arg_parser.parse_args()

# Directories to store data files
ROOT_DIR = C["ROOT_DIR"]
UNIPROT_DIR = os.path.join(ROOT_DIR, "external/UNIPROT")
PDB_DIR = os.path.join(ROOT_DIR, "external/PDB")
# CATH_DIR = os.path.join(ROOT_DIR, "external/mappings/CATH")
# SIFTS_DIR = os.path.join(ROOT_DIR, "external/mappings/SIFTS")
# GO_DIR = os.path.join(ROOT_DIR, "external/mappings/GO")
# CIF_DIR = os.path.join(ROOT_DIR, "external/CIFFILES")
# MAP_DIR = os.path.join(ROOT_DIR, "external/mappings")


# Get list of valid PDBids
print("Getting list of PDB ids")
PDBIDS = [pdbid.strip().lower() for pdbid in open(ARGS.pdb_list_file)]
print(PDBIDS)

# Perform GraphQL query
graphql_url = "https://data.rcsb.org/graphql"
query = open(os.path.join(ROOT_DIR, "external/graphql_query.tpl")).read()
query = Template(query)
QUERY = {
    "query": query.substitute(entry_ids_array=','.join(['"%s"' %s for s in PDBIDS]))
}
print("Performing GraphQL request at {}".format(graphql_url))
req = requests.post(graphql_url, json=QUERY)
print(req.status_code)
ENTRY_DATA = req.json()

POLYMER_ENTITY_INSTANCE_DATA = {}
# Loop over every entry returned by query
for entry in ENTRY_DATA["data"]["entries"]:
    pdb_id = entry["rcsb_id"]
    
    # loop over polymer entities
    for poly_entity in entry["polymer_entities"]:
        if poly_entity["entity_poly"]["rcsb_entity_polymer_type"] != "Protein":
            continue
        
        # get annotations
        polymer_id = poly_entity["rcsb_id"]
        clusters = {c["identity"]: c["cluster_id"] for c in poly_entity["rcsb_cluster_membership"]}
        uniprot_accessions = [up["rcsb_id"] for up in poly_entity["uniprots"]]
        
        # loop over polymer instances
        for poly_entity_instance in poly_entity["polymer_entity_instances"]:
            polymer_instance_id = poly_entity_instance["rcsb_id"]
            asym_id = poly_entity_instance["rcsb_polymer_entity_instance_container_identifiers"]["asym_id"] # PDB assigned chain id
            auth_asym_id = poly_entity_instance["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"] # author assigned chain id
            
            # create chain features dictionary
            POLYMER_ENTITY_INSTANCE_DATA[polymer_instance_id] = {
                "cath": {
                    "H": [],
                    "T": [],
                    "A": [],
                    "C": [],
                },
                "uniprot": {
                    "accession": uniprot_accessions,
                    "names": [],
                    "organism": 'N/A',
                    "keywords":{
                        "kw_id": [],
                        "category": [],
                        "name": []
                    }
                },
                "go": {
                    "molecular_function": [],
                    "biological_process": [],
                    "cellular_component": [],
                },
                "clusters": clusters,
                "chain_id": auth_asym_id,
                "asym_id": asym_id
            }
            
            # add instance features
            for annotation in poly_entity_instance["rcsb_polymer_instance_annotation"]:
                if annotation["type"] == "CATH":
                    cath = annotation["annotation_id"].split('.')
                    POLYMER_ENTITY_INSTANCE_DATA[polymer_instance_id]["cath"]["C"].append(cath[0])
                    POLYMER_ENTITY_INSTANCE_DATA[polymer_instance_id]["cath"]["A"].append(cath[1])
                    POLYMER_ENTITY_INSTANCE_DATA[polymer_instance_id]["cath"]["T"].append(cath[2])
                    POLYMER_ENTITY_INSTANCE_DATA[polymer_instance_id]["cath"]["H"].append(cath[3])

# Retrieve UniProt Data
if ARGS.uniprot_file:
    with open(ARGS.uniprot_file) as FH:
        UNIPROT_DATA_BACKUP = json.load(FH)
else:
    UNIPROT_DATA_BACKUP = None
UNIPROT_DATA = {}
uniprot_url = "https://rest.uniprot.org/uniprotkb/{}.json"
GO_categories = {
    "F": "molecular_function",
    "C": "cellular_component",
    "P": "biological_process"
}
for polymer_instance_id in POLYMER_ENTITY_INSTANCE_DATA:
    polymer_instance = POLYMER_ENTITY_INSTANCE_DATA[polymer_instance_id]
    accessions = polymer_instance["uniprot"]["accession"]
    
    for accession in accessions:
        # download if not already
        if accession in UNIPROT_DATA:
            up = UNIPROT_DATA[accession]
        else:
            req = requests.get(uniprot_url.format(accession))
            if req.status_code == "200":
                up = req.json()
            else:
                # fallback to backup if given
                if UNIPROT_DATA_BACKUP and accession in UNIPROT_DATA_BACKUP:
                    up = UNIPROT_DATA_BACKUP[accession]
                else:
                    continue # no source could be found for this accesssion
            UNIPROT_DATA[accession] = up
        
        # add uniprot info
        polymer_instance["uniprot"]["organism"] = up["organism"]["scientificName"]
        polymer_instance["uniprot"]["names"].append(up["proteinDescription"]["recommendedName"]["fullName"]["value"])
        if "alternativeNames" in up["proteinDescription"]:
            for name in up["proteinDescription"]["alternativeNames"]:
                polymer_instance["uniprot"]["names"].append(name["fullName"]["value"])
        for kw in up["keywords"]:
            polymer_instance["uniprot"]["keywords"]["kw_id"].append(kw["id"])
            polymer_instance["uniprot"]["keywords"]["category"].append(kw["category"])
            polymer_instance["uniprot"]["keywords"]["name"].append(kw["name"])
        for cr in up["uniProtKBCrossReferences"]:
            if cr["database"] == "GO":
                go_id = cr["id"]
                cat, description = cr["properties"][0]["value"].split(':')
                polymer_instance["go"][GO_categories[cat]].append({"description":description, "GO_ID":go_id})

print(UNIPROT_DATA)
print(UNIPROT_DATA_BACKUP)
# Save UniProt data to disk
if ARGS.refresh_uniprot:
    if UNIPROT_DATA_BACKUP is None:
        UNIPROT_DATA_BACKUP = {}
    UNIPROT_DATA_BACKUP.update(UNIPROT_DATA)
    print(UNIPROT_DATA_BACKUP)
    with open(os.path.join(UNIPROT_DIR, "uniprot_data.json"), "w") as FH:
        FH.write(json.dumps(UNIPROT_DATA_BACKUP))

# Write cluster mappings to file
cluster_url = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{}.txt"
CLUSTERS = ["30", "40", "50", "70", "90", "95", "100"]
if ARGS.download_clusters:
    for cluster in CLUSTERS:
        # download cluster file
        req = requests.get(cluster_url.format(cluster))
        clusters = req.text
        cluster_map = {}
        for i, cl in enumerate(clusters.splitlines()):
            cnum = i + 1
            for poly_ent in cl.strip().split():
                if poly_ent[0:2] == "AF":
                    continue
                if poly_ent[0:2] == "MA":
                    continue
                cluster_map[poly_ent] = cnum
        
        # save map
        path = os.path.join(PDB_DIR, "cluster-{}.pkl".format(cluster))
        with open(path, "wb") as FH:
            pickle.dump(cluster_map, FH)
exit(0)

# write PDBID data to file
print("Writing PDB id info to file")
for pdbid in PDBIDS:
    d = pdbid[-1]
    for chain in PDBIDS[pdbid]:
        ckey = "{}_{}".format(pdbid.upper(), chain)
        # add sequence clusters
        for cluster in CLUSTERS:
            #path = os.path.join(PDB_DIR, d, "{}.{}_{}.xml".format(pdbid, chain, cluster))
            #if(os.access(path, os.R_OK)):
            #    REP = open(path)
            #    data = xmltodict.parse(REP.read())
            #    REP.close()
            #    if(data['representatives']):
            #        PDBIDS[pdbid][chain]["clusters"][cluster] = data['representatives']['pdbChain']['@name']
            #    else:
            #        PDBIDS[pdbid][chain]["clusters"][cluster] = 'N/A'
            #else:
            #    PDBIDS[pdbid][cid]["clusters"][cluster] = 'N/A'
            if(ckey in CLUSTER_MAP[cluster]):
                PDBIDS[pdbid][chain]["clusters"][cluster] = CLUSTER_MAP[cluster][ckey]
            else:
                PDBIDS[pdbid][chain]["clusters"][cluster] = 'N/A'
        
        # check for empty CATH data
        if(len(PDBIDS[pdbid][chain]["cath"]["H"]) == 0):
            PDBIDS[pdbid][chain]["cath"]["H"].append('N/A')
            PDBIDS[pdbid][chain]["cath"]["T"].append('N/A')
            PDBIDS[pdbid][chain]["cath"]["A"].append('N/A')
            PDBIDS[pdbid][chain]["cath"]["C"].append('N/A')
        
        # check for empty GO data
        if(
            len(PDBIDS[pdbid][chain]["go"]["molecular_function"]) == 0 and
            len(PDBIDS[pdbid][chain]["go"]["biological_process"]) == 0 and
            len(PDBIDS[pdbid][chain]["go"]["cellular_component"]) == 0
        ):
            PDBIDS[pdbid][chain]["go"]["molecular_function"].append({"description": 'N/A', "GO_ID": 'N/A'})
            PDBIDS[pdbid][chain]["go"]["biological_process"].append({"description": 'N/A', "GO_ID": 'N/A'})
            PDBIDS[pdbid][chain]["go"]["cellular_component"].append({"description": 'N/A', "GO_ID": 'N/A'})
        
        # check for empty Uniprot Data
        if(len(PDBIDS[pdbid][chain]['uniprot']['accession']) == 0):
            PDBIDS[pdbid][chain]['uniprot']['accession'].append('N/A')
        
        # remove sets
        del PDBIDS[pdbid][chain]["go"]["seen"]
        del PDBIDS[pdbid][chain]["cath"]["seen"]
        del PDBIDS[pdbid][chain]["uniprot"]["seen"]
    FH = open(os.path.join(MAP_DIR, "{}/{}.json".format(d, pdbid)), "w")
    FH.write(json.dumps(PDBIDS[pdbid]))
    FH.close()

