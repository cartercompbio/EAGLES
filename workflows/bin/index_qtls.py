#!/usr/bin/env python

import pickle
import argparse       

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get gene-associated eQTLs")
    parser.add_argument("--input", required=True, help="Path to qtl table")
    parser.add_argument("--output", required=True, help="Path to output qtl table pickle file")

    args = parser.parse_args()

    index_gtex_file(args.input, args.output)