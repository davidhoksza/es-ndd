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
    for f_name in glob.glob(dir_name + "/*.txt"):
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
        for unp_id in genes[gene_name]:
            genes[gene_name][unp_id] = list(genes[gene_name][unp_id])

    return genes


def reorganize(data):
    genes = set()
    for type in DS_TYPE:
        if type in data:
            genes = genes.union(data[type])

    data_r = {}
    for gene in genes:
        data_r[gene] = {}
        for type in DS_TYPE:
            if type not in data:
                continue
            if gene in data[type]:
                data_r[gene][type.name] = data[type][gene]

    return data_r


def main():
    data = {}
    for type in DS_TYPE:
        data[type] = process_dir_single_structure(type)
    print(data)
    with open("list.json", "w") as f:
        f.write(json.dumps(reorganize(data), indent=2))

if __name__ == "__main__":
    main()

