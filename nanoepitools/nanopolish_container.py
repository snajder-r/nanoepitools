import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from typing import Union, List

from nanoepitools.nanopolish_calls import SparseMethylationMatrixContainer


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


class MethlyationValuesContainer:
    
    def __init__(self, chromosome_container, start, end):
        self.chromosome = chromosome_container
        self.start = start
        self.end = end
    
    def get_read_names_unique(self):
        read_name_ds = self.chromosome.h5group['read_name'][self.start:self.end]
        reads_names, idx = np.unique(read_name_ds, return_index=True)
        reads_names = reads_names[np.argsort(idx)]
        return reads_names
    
    def get_ranges_unique(self):
        ranges_ds = self.chromosome.h5group['range'][self.start:self.end]
        # Since they are already pre-sorted, this is fast way to get uniques:
        diff = np.ones_like(ranges_ds)
        diff[1:] = ranges_ds[1:] - ranges_ds[:-1]
        idx = diff[:, 0].nonzero()[0] & diff[:, 1].nonzero()[0]
        return ranges_ds[idx, :]
    
    def get_ranges(self):
        return self.chromosome.h5group['range'][self.start:self.end, :]
    
    def get_llrs(self):
        return self.chromosome.h5group['llr'][self.start:self.end]
    
    def get_read_names(self):
        return self.chromosome.h5group['read_name'][self.start:self.end]
    
    def to_sparse_methylation_matrix(self) -> SparseMethylationMatrixContainer:
        read_names = [r.decode() for r in self.get_read_names_unique()]
        genomic_ranges = self.get_ranges_unique()
        
        coord_to_index_dict = {genomic_ranges[i, 0]: i for i in
                               range(len(genomic_ranges))}
        
        met_matrix = np.zeros((len(read_names), len(genomic_ranges)))
        read_dict = {read_names[i]: i for i in range(len(read_names))}
        
        range_ds = self.get_ranges()
        read_name_ds = self.get_read_names()
        llr_ds = self.get_llrs()
        for i in range(len(range_ds)):
            read_i = read_dict[read_name_ds[i].decode()]
            range_i = coord_to_index_dict[range_ds[i, 0]]
            met_matrix[read_i, range_i] = llr_ds[i]
        met_matrix = sp.csc_matrix(met_matrix)
        return SparseMethylationMatrixContainer(met_matrix, read_names,
                                                genomic_ranges[:, 0],
                                                genomic_ranges[:, 1])


class ChromosomeContainer:
    
    def __init__(self, chromosome_group: h5py.Group, chunk_size: int):
        self.h5group = chromosome_group
        self.chunk_size = chunk_size
    
    def __len__(self):
        return len(self.h5group['range'])
    
    def get_number_of_chunks(self):
        num_chunks = len(self) // self.chunk_size
        if len(self) % self.chunk_size != 0:
            num_chunks += 1
        return num_chunks
    
    def get_chunk_ids(self):
        return [i for i in range(self.get_number_of_chunks())]
    
    def _seek_overlap_ranges_backwards(self, chunk_id, start_value=-1):
        if start_value == -1:
            start_value = self.h5group['range'][self.chunk_size * chunk_id, 0]
        
        starts = self.h5group['range'][(self.chunk_size * chunk_id):(
                self.chunk_size * (chunk_id + 1)), 0]
        matches = np.arange(self.chunk_size)[starts == start_value]
        
        if len(matches) == 0:
            # Nothing in this chunk, return beginning of the chunk we came from
            return self.chunk_size * (chunk_id + 1)
        
        if matches[0] == 0 and chunk_id > 0:
            # All of this chunk is the same range, we need to go deeper
            return self._seek_overlap_ranges_backwards(chunk_id - 1,
                                                       start_value=start_value)
        
        # Part of this chunk has entries for this start position
        return self.chunk_size * chunk_id + matches[0]
    
    def _seek_overlap_ranges_forwards(self, chunk_id, end_value=-1):
        last = min(len(self), self.chunk_size * (chunk_id + 1) - 1)
        
        if end_value == -1:
            end_value = self.h5group['range'][last, 0]
        
        ends = self.h5group['range'][
               (self.chunk_size * chunk_id):(last+1),
               0]
        
        matches = np.arange(len(ends))[ends == end_value]
        
        if len(matches) == 0:
            # Nothing in this chunk, return end of the chunk we came from
            return self.chunk_size * chunk_id - 1
        
        if matches[
            -1] == self.chunk_size - 1 and chunk_id < \
                self.get_number_of_chunks() - 1:
            # All of this chunk is the same range, we need to go deeper
            return self._seek_overlap_ranges_forwards(chunk_id + 1,
                                                      end_value=end_value)
        
        # Part of this chunk has entries for this end position
        return self.chunk_size * chunk_id + matches[-1]
    
    def get_chunk(self, chunk_id: int, overlap=True):
        """
        Returns the values of said chunk, and, if overlap=True, includes values
        of neighboring chunks if they are in the same genomic ranges, such as
        to avoid having a subset of reads of one location in one chunk and the
        rest in the other
        :param chunk_id: The chunk id (see get_chunk_ids)
        :param overlap: Whether to look for same-region locations in
        neighboring chunks
        :return: MethlyationValuesContainer
        """
        if overlap:
            earliest_pos = self._seek_overlap_ranges_backwards(chunk_id)
            latest_pos = self._seek_overlap_ranges_forwards(chunk_id) + 1
        else:
            earliest_pos = self.chunk_size * chunk_id
            latest_pos = min(self.chunk_size * (chunk_id + 1), len(self))
        
        return MethlyationValuesContainer(self, earliest_pos, latest_pos)


class MetcallH5Container:
    
    def __init__(self, h5filepath: Union[str, Path], mode: str = 'r',
                 chunk_size=100000):
        self.h5filepath = h5filepath
        self.mode = mode
        self.chunk_size = chunk_size
        self.h5_fp = None
        self.chrom_container_cache = {}
    
    def __enter__(self):
        self.h5_fp = h5py.File(self.h5filepath, mode=self.mode)
        return self
    
    def __exit__(self, exittype, exitvalue, traceback):
        self.h5_fp.close()
    
    def add_to_h5_file(self, cur_df):
        main_group = self.h5_fp.require_group('chromosomes')
        
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
    
    def get_chromosomes(self) -> List[str]:
        return [str(k) for k in self.h5_fp['chromosomes'].keys()]
    
    def __getitem__(self, chromosome: str):
        if chromosome not in self.h5_fp['chromosomes'].keys():
            return ValueError('No data for this chromosome in container')
        if chromosome in self.chrom_container_cache.keys():
            return self.chrom_container_cache[chromosome]
        else:
            ret = ChromosomeContainer(self.h5_fp['chromosomes'][chromosome],
                                      self.chunk_size)
            self.chrom_container_cache[chromosome] = ret
            return ret
