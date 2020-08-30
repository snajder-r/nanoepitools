from pathlib import Path
from typing import Union

import h5py
import numpy as np
import pandas as pd


def create_or_extend(parent_group, name, shape, data, **kwargs):
    if name not in parent_group.keys():
        parent_group.create_dataset(name=name, shape=shape, data=data, **kwargs)
    else:
        num_data = data.shape[0] if hasattr(data, 'shape') else len(data)
        ds = parent_group[name]
        old_shape = ds.shape
        new_shape = (old_shape[i] + (num_data if i == 0 else 0) for i in
                     range(len(old_shape)))
        ds.resize(new_shape)
        
        ds[old_shape[0]:] = data
        
        print("Extended from ", old_shape, " to ", ds.shape)


def argsort_ranges(chrom_group):
    return np.argsort(chrom_group['range'][:, 0], kind='mergesort')

    
class MetcallH5Container:
    
    def __init__(self, outfile: Union[str, Path], mode: str = 'r',
                 chunk_size=100000):
        self.outfile = outfile
        self.mode = mode
        self.chunk_size = chunk_size
        self.out_fp = None
    
    def __enter__(self):
        self.out_fp = h5py.File(self.outfile, mode=self.mode)
        return self
    
    def __exit__(self, exittype, exitvalue, traceback):
        self.out_fp.close()
    
    def add_to_h5_file(self, cur_df):
        main_group = self.out_fp.require_group('chromosomes')
        
        for chrom in set(cur_df['chromosome']):
            print(chrom)
            
            chrom_calls = cur_df.loc[cur_df['chromosome'] == chrom]
            n = chrom_calls.shape[0]
            read_names = [read.encode() for read in chrom_calls['read_name']]
            read_name_len = len(read_names[0])
            assert all([len(read) for read in read_names])
            
            chrom_group = main_group.require_group(chrom)
            
            create_or_extend(parent_group=chrom_group, name='range',
                             shape=(n, 2), dtype=int,
                             data=chrom_calls[['start', 'end']],
                             compression='gzip', chunks=(self.chunk_size, 2),
                             maxshape=(None, 2))
            create_or_extend(parent_group=chrom_group, name='llr', shape=(n,),
                             dtype=float, data=chrom_calls['log_lik_ratio'],
                             compression='gzip', chunks=(self.chunk_size,),
                             maxshape=(None,))
            create_or_extend(parent_group=chrom_group, name='read_name',
                             shape=(n,), dtype='S%d' % read_name_len,
                             data=read_names, compression='gzip',
                             chunks=(self.chunk_size,), maxshape=(None,))
            
            sort_order = argsort_ranges(chrom_group)
            chrom_group['range'][:] = np.array(chrom_group['range'])[sort_order]
            chrom_group['llr'][:] = np.array(chrom_group['llr'])[sort_order]
            chrom_group['read_name'][:] = np.array(chrom_group['read_name'])[
                sort_order]
    
    def parse_and_add_nanopolish_file(self, nanopolish_file: Union[str, Path]):
        cur_df = pd.read_csv(nanopolish_file, sep='\t',
                             dtype={'chromosome': str})
        self.add_to_h5_file(cur_df)
