import glob
import os
import enum
import json

class DS_TYPE(enum.Enum):
    PDB = 1    
    ALPHAFOLD = 2
    PDB_MULTI = 3


def get_dir_name(type: DS_TYPE):
    if type == DS_TYPE.PDB:
        return "pdb"
    if type == DS_TYPE.PDB_MULTI:
        return "pdb-multi"    
    if type == DS_TYPE.ALPHAFOLD:
        return "alphafold"


def add_to_genes_list(gl, gene_name, unp_id, pdb_id):
    if gene_name not in gl:
        gl[gene_name] = {}

    if unp_id not in gl[gene_name]:
        gl[gene_name][unp_id] = set()

    gl[gene_name][unp_id].add(pdb_id)


def process_dir_single_structure(type: DS_TYPE):

    genes = {}    
    dir_name = get_dir_name(type)
    
    preferred_structures = {} # dictionary of gene and pdb id to be shown first in the list of available structures in the web app
    if type == DS_TYPE.PDB or type == DS_TYPE.PDB_MULTI:
        with open("{}/preferred_starting_pdbids.txt".format(dir_name)) as f:
            f.readline()
            for line in f.readlines():
                line = line.strip()
                if line == "":
                    continue
                gene, pdb_id = line.split('\t')
                preferred_structures[gene] = pdb_id        
    
    
    for f_name in glob.glob(dir_name + "/*.txt"):
        if 'preferred_starting' in f_name: 
            # preferred_starting_pdbids_multimer_dataset.txt
            continue
        if type == DS_TYPE.PDB:
            aux = os.path.basename(f_name).replace('Experimentally_solved_sinlge_', '').replace('.txt', '')
            pos = aux.find('_')
            gene_name = aux[:pos]
            pdb_id = aux[pos+1:]
            unp_ids = set()        
        elif type == DS_TYPE.ALPHAFOLD:
            aux = os.path.basename(f_name).replace('Alpha_fold_', '').replace('.txt', '')
            gene_name = aux
            pdb_id = aux
            unp_ids = set()
        elif type == DS_TYPE.PDB_MULTI:
            aux = os.path.basename(f_name).replace('Experimentally_solved_complex_', '').replace('.txt', '')
            pdb_id = aux

        with open(f_name) as f:
            f.readline()
            for line in f.readlines():
                if type == DS_TYPE.PDB_MULTI:
                    sLine = line.split('\t')
                    unp_id = sLine[-3].strip()
                    gene_name = sLine[-2].strip(" \"")
                    chain = sLine[-1].strip()
                    add_to_genes_list(gl=genes, gene_name=gene_name, unp_id=unp_id, pdb_id=pdb_id+"_"+chain)
                else:
                    unp_id = line.split('\t')[-1].strip()
                    unp_ids.add(unp_id)

        if type != DS_TYPE.PDB_MULTI:
            assert len(unp_ids) == 1
            unp_id = unp_ids.pop()
            add_to_genes_list(gl=genes, gene_name=gene_name, unp_id=unp_id, pdb_id=pdb_id)

    for gene_name in genes:
        preferred_pdb_id = preferred_structures[gene_name] if gene_name in preferred_structures else None        
        if type == DS_TYPE.PDB and preferred_pdb_id: 
            #for non PDB (not PDB_multi) the id is 2fo0A and needs to be transformed to 2fo0_A
            if '_' not in preferred_pdb_id:
                preferred_pdb_id = preferred_pdb_id[:-1] + '_' + preferred_pdb_id[-1]                
            
        for unp_id in genes[gene_name]:
            # sort while starting with preferred pdb_id
            pdb_list = sorted(list(genes[gene_name][unp_id]))
            if preferred_pdb_id:
                preferred_pos = 0
                for i in range(len(pdb_list)):
                    if pdb_list[i].startswith(preferred_pdb_id):
                        preferred_pos = i
                        break
                pdb_list.insert(0, pdb_list.pop(preferred_pos))
            genes[gene_name][unp_id] = pdb_list

    return genes


def reorganize(data):
    f = open("ndd_genes.txt")
    genes_filter = [l.strip().strip('"') for l in f.readlines()]
        
    genes = set()
    for type in DS_TYPE:
        if type in data:
            genes = genes.union(data[type])

    data_r = {}
    for gene in genes:
        if gene not in genes_filter:
            continue
        data_r[gene] = {}
        for type in DS_TYPE:
            if type not in data:
                continue
            if gene in data[type]:
                data_r[gene][type.name] = data[type][gene]
                
    print('Not used genes from the filter: ', list(set(genes_filter)-set(data_r.keys())))

    return data_r

   
def main():
    data = {}
    for type in DS_TYPE:
        data[type] = process_dir_single_structure(type)    
    # print(data)
    with open("list.json", "w") as f:
        f.write(json.dumps(reorganize(data), indent=2))

if __name__ == "__main__":    
    main()

