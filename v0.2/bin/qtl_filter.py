#!/usr/bin/env python

import pickle
import os
import argparse
from functools import reduce

def index_gtex_file(filename, index_file):
    
    if filename.endswith('.csv'):
        sep = ','
    elif filename.endswith('.tsv'):
        sep = '\t'
    else:
        assert False
    
    id_ranges = {}
    
    with open(filename, 'rb') as f:
        header_end = f.tell() + len(f.readline())
        
        current_id = None
        start_byte = header_end
        
        while True:
            line_start = f.tell()
            line = f.readline()
            if not line:
                break
            
            row_id = line.decode().split(sep)[0]
            
            if row_id != current_id:
                if current_id is not None:
                    try:
                        id_ranges[current_id]
                        assert False
                    except KeyError:
                        id_ranges[current_id] = (start_byte, line_start)
                    except AssertionError:
                        raise ValueError('file improperly formatted, file must be sorted by identifiers in column 1')
                current_id = row_id
                start_byte = line_start
        
        if current_id is not None:
            id_ranges[current_id] = (start_byte, f.tell())
    
    with open(index_file, 'wb') as idx:
        pickle.dump(id_ranges, idx)
        
def load_by_id(filename, index_file, gene, col_name = 'SNP'):
    if filename.endswith('.csv'):
        sep = ','
    elif filename.endswith('.tsv'):
        sep = '\t'
    else:
        raise ValueError(f'filename must be .csv or .tsv, not {filename}')
        
    with open(index_file, 'rb') as idx:
        id_ranges = pickle.load(idx)
        
    if gene not in id_ranges:
        return []
        
    start_byte,end_byte = id_ranges[gene]

    with open(filename, 'rb') as infile:
        l = next(infile).decode().strip().split(sep)
        row_pos = l.index(col_name)

        infile.seek(start_byte)

        values = []
        while infile.tell() < end_byte:
            line = infile.readline().decode().strip()
            if line:
                values.append(line.split(sep)[row_pos])

        return values
    
def load_by_tis(gene, outfile, **kwargs):
    # ex tissues: ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Cells_Cultured_fibroblasts', 'Cells_EBV-transformed_lymphocytes', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Kidney_Cortex', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood']
    try:
        tis_list = kwargs['tis_list']
    except:
        tis_list = ['Colon_Sigmoid', 'Colon_Transverse']
        
    try:
        merge_type = kwargs['merge_type']
    except:
        merge_type = 'any'
        
    try:
        qtl_folder = kwargs['qtl_folder']
    except:
        raise ValueError('qtl_folder kwarg missing')
        
    try:
        index_folder = kwargs['index_folder']
    except:
        raise ValueError('index_folder kwarg missing')
        
    try:
        col_name = kwargs['col_name']
    except:
        col_name = 'SNP'
        
        
    if merge_type == 'any':
        func = lambda a,b: a|b
    elif merge_type == 'all':
        func = lambda a,b: a&b
        
    snp_sets = [set(load_by_id(f'{qtl_folder}/{t}.tsv', f'{index_folder}/{t}.pkl', gene)) for t in tis_list]
    merged_set = reduce(func, snp_sets)
    
    if len(merged_set)>0:
        with open(outfile, 'w') as outfile:
            outfile.write('\n'.join(merged_set))        
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get gene-associated eQTLs")
    parser.add_argument("--qtl-folder", required=True, help="Path to qtl table")
    parser.add_argument("--index-folder", required=True, help="Path to qtl table index file")
    parser.add_argument("--gene", required=True, help="gene identifier for desired eQTLs")
    parser.add_argument("--tis", nargs='+', required=True, help="list of tissues to consider")
    parser.add_argument("--output", required=True, help="Path to output eqtl list")

    args = parser.parse_args()

    # Load data
    load_by_tis(args.gene, args.output, qtl_folder = args.qtl_folder, index_folder = args.index_folder, tis_list = args.tis)
