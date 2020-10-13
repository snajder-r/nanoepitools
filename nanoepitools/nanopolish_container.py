from __future__ import annotations
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import logging
from pathlib import Path
from typing import Union, List, Tuple, Dict

from nanoepitools.nanopolish_calls import SparseMethylationMatrixContainer


def unique_genomic_range(genomic_ranges: nd.array):
    """
    :param genomic_ranges: List of tuples, must be PRESORTED
    :return:
    """
    diff = np.ones_like(genomic_ranges)
    diff[1:] = genomic_ranges[1:] - genomic_ranges[:-1]
    idx = diff[:, 0].nonzero()[0] & diff[:, 1].nonzero()[0]
    return genomic_ranges[idx, :]


def argsort_ranges(chrom_group):
    return np.argsort(chrom_group["range"][:, 0], kind="mergesort")


def create_sparse_matrix_from_samples(
    sample_llrs: Dict[str:MethlyationValuesContainer],
) -> SparseMethylationMatrixContainer:
    samples = list(sample_llrs.keys())
    read_names = {
        s: [r.decode() for r in sample_llrs[s].get_read_names_unique()] for s in samples
    }
    genomic_ranges = {
        s: [r for r in sample_llrs[s].get_ranges_unique()] for s in samples
    }

    sample_assignment = [s for s in samples for _ in read_names[s]]
    read_names = [r for s in samples for r in read_names[s]]
    genomic_ranges = [r for s in samples for r in genomic_ranges[s]]
    genomic_ranges = np.array(sorted(genomic_ranges, key=lambda x: x[0] * 10e10 + x[1]))
    genomic_ranges = unique_genomic_range(genomic_ranges)

    coord_to_index_dict = {genomic_ranges[i, 0]: i for i in range(len(genomic_ranges))}

    met_matrix = np.zeros((len(read_names), len(genomic_ranges)))
    read_dict = {read_names[i]: i for i in range(len(read_names))}

    for sample, llrs in sample_llrs.items():
        range_ds = llrs.get_ranges()
        read_name_ds = llrs.get_read_names()
        llr_ds = llrs.get_llrs()
        for i in range(len(range_ds)):
            read_i = read_dict[read_name_ds[i].decode()]
            range_i = coord_to_index_dict[range_ds[i, 0]]
            met_matrix[read_i, range_i] = llr_ds[i]

    met_matrix = sp.csc_matrix(met_matrix)
    return SparseMethylationMatrixContainer(
        met_matrix,
        read_names,
        genomic_ranges[:, 0],
        genomic_ranges[:, 1],
        read_samples=sample_assignment,
    )


class MethlyationValuesContainer:
    def __init__(self, chromosome_container, start, end):
        self.chromosome = chromosome_container
        self.start = start
        self.end = end

    def get_read_names_unique(self):
        read_name_ds = self.chromosome.h5group["read_name"][self.start : self.end]
        reads_names, idx = np.unique(read_name_ds, return_index=True)
        reads_names = reads_names[np.argsort(idx)]
        return reads_names

    def get_ranges_unique(self):
        ranges_ds = self.chromosome.h5group["range"][self.start : self.end]
        return unique_genomic_range(ranges_ds)

    def get_ranges(self):
        return self.chromosome.h5group["range"][self.start : self.end, :]

    def get_llrs(self):
        return self.chromosome.h5group["llr"][self.start : self.end]

    def get_read_names(self):
        return self.chromosome.h5group["read_name"][self.start : self.end]

    def get_read_groups(self, group_key):
        return self.chromosome.h5group["read_groups"][group_key][self.start : self.end]

    def to_sparse_methylation_matrix(
        self, read_groups_key: str = None
    ) -> SparseMethylationMatrixContainer:
        read_names = [r.decode() for r in self.get_read_names_unique()]
        genomic_ranges = self.get_ranges_unique()

        coord_to_index_dict = {
            genomic_ranges[i, 0]: i for i in range(len(genomic_ranges))
        }

        met_matrix = np.zeros((len(read_names), len(genomic_ranges)))
        read_dict = {read_names[i]: i for i in range(len(read_names))}

        range_ds = self.get_ranges()
        read_name_ds = self.get_read_names()
        llr_ds = self.get_llrs()

        if read_groups_key is not None:
            read_groups_ds = self.get_read_groups(read_groups_key)
            read_samples_dict = {
                read_name_ds[i].decode(): read_groups_ds[i]
                for i in range(len(read_groups_ds))
            }
            read_samples = [read_samples_dict[r] for r in read_names]
        else:
            read_samples = None

        for i in range(len(range_ds)):
            read_i = read_dict[read_name_ds[i].decode()]
            range_i = coord_to_index_dict[range_ds[i, 0]]
            met_matrix[read_i, range_i] = llr_ds[i]
        met_matrix = sp.csc_matrix(met_matrix)
        return SparseMethylationMatrixContainer(
            met_matrix,
            read_names,
            genomic_ranges[:, 0],
            genomic_ranges[:, 1],
            read_samples=read_samples,
        )


class ChromosomeContainer:
    def __init__(self, chromosome_group: h5py.Group, chunk_size: int):
        self.h5group = chromosome_group
        self.chunk_size = chunk_size

    def __len__(self):
        return len(self.h5group["range"])

    def get_number_of_chunks(self):
        num_chunks = len(self) // self.chunk_size
        if len(self) % self.chunk_size != 0:
            num_chunks += 1
        return num_chunks

    def get_chunk_ids(self):
        return [i for i in range(self.get_number_of_chunks())]

    def _seek_overlap_ranges_backwards(self, chunk_id, start_value=-1):
        last = min(len(self), self.chunk_size * (chunk_id + 1)) - 1
        if start_value == -1:
            start_value = self.h5group["range"][self.chunk_size * chunk_id, 0]

        starts = self.h5group["range"][(self.chunk_size * chunk_id) : last, 0]
        matches = np.arange(len(starts))[starts == start_value]

        if len(matches) == 0:
            # Nothing in this chunk, return beginning of the chunk we came from
            return self.chunk_size * (chunk_id + 1)

        if matches[0] == 0 and chunk_id > 0:
            # All of this chunk is the same range, we need to go deeper
            return self._seek_overlap_ranges_backwards(
                chunk_id - 1, start_value=start_value
            )

        # Part of this chunk has entries for this start position
        return self.chunk_size * chunk_id + matches[0]

    def _seek_overlap_ranges_forwards(self, chunk_id, end_value=-1):
        last = min(len(self), self.chunk_size * (chunk_id + 1)) - 1

        if end_value == -1:
            end_value = self.h5group["range"][last, 0]

        ends = self.h5group["range"][(self.chunk_size * chunk_id) : (last + 1), 0]

        matches = np.arange(len(ends))[ends == end_value]

        if len(matches) == 0:
            # Nothing in this chunk, return end of the chunk we came from
            return self.chunk_size * chunk_id - 1

        if (
            matches[-1] == self.chunk_size - 1
            and chunk_id < self.get_number_of_chunks() - 1
        ):
            # All of this chunk is the same range, we need to go deeper
            return self._seek_overlap_ranges_forwards(chunk_id + 1, end_value=end_value)

        # Part of this chunk has entries for this end position
        return self.chunk_size * chunk_id + matches[-1]

    def get_chunk(self, chunk_id: int, overlap=True):
        """Returns the values of said chunk, and, if overlap=True,
        includes values of neighboring chunks if they are in the same
        genomic ranges, such as to avoid having a subset of reads of one
        location in one chunk and the rest in the other.

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

    def create_chunk_index(self, force_update=False):
        if "chunk_ranges" in self.h5group.keys() and not force_update:
            return

        index = np.zeros((self.get_number_of_chunks(), 2))
        num_ranges = self.h5group["range"].shape[0]
        for chunk_i, start_i in enumerate(range(0, num_ranges, self.chunk_size)):
            end_i = min(num_ranges - 1, start_i + self.chunk_size)
            index[chunk_i, 0] = self.h5group["range"][start_i, 0]
            index[chunk_i, 1] = self.h5group["range"][end_i, 1]

        if "chunk_ranges" in self.h5group.keys():
            self.h5group["chunk_ranges"].resize(index.shape)
            self.h5group["chunk_ranges"][:] = index
        else:
            self.h5group.create_dataset(
                name="chunk_ranges", data=index, dtype=int, maxshape=(None, 2)
            )
        self.h5group.attrs["chunk_size"] = self.chunk_size

    def get_values_in_range(self, genomic_start, genomic_end):
        if "chunk_size" not in self.h5group.attrs.keys():
            raise ValueError(
                "Random access to ranges only allowed if index exists. Call create_chunk_index"
            )
        index_chunk_size = self.h5group.attrs["chunk_size"]
        index = self.h5group["chunk_ranges"][:]

        # First find the right chunk for start and end

        chunk_indices = np.arange(len(index))[
            (index[:, 0] < genomic_end) & (genomic_start <= index[:, 1])
        ]

        if len(chunk_indices) == 0:
            # If no chunk contains these values
            return None

        start_chunk = chunk_indices[0]
        end_chunk = chunk_indices[-1]

        # Find precise start point
        start_index = start_chunk * index_chunk_size
        start_chunk_start = start_chunk * index_chunk_size
        start_chunk_end = min(len(self) - 1, (start_chunk + 1) * index_chunk_size)
        start_chunk_ranges = self.h5group["range"][start_chunk_start:start_chunk_end, :]
        start_in_range_indices = np.arange(len(start_chunk_ranges))[
            start_chunk_ranges[:, 1] >= genomic_start
        ]
        if len(start_in_range_indices) > 0:
            # Add index of first value that is in the range
            start_index += start_in_range_indices[0]

        # Find precise end point
        end_index = end_chunk * index_chunk_size
        end_chunk_start = end_chunk * index_chunk_size
        end_chunk_end = min(len(self) - 1, (end_chunk + 1) * index_chunk_size)
        end_chunk_ranges = self.h5group["range"][end_chunk_start:end_chunk_end, :]
        end_oor_indices = np.arange(len(end_chunk_ranges))[
            end_chunk_ranges[:, 0] >= genomic_end
        ]
        if len(end_oor_indices) > 0:
            # Add index of first value that is out of range
            end_index += end_oor_indices[0]
        else:
            # If all values in the chunk are in the range
            end_index = min(len(self), end_index + index_chunk_size)

        return MethlyationValuesContainer(self, start_index, end_index)


class MetcallH5Container:
    def __init__(
        self, h5filepath: Union[str, Path], mode: str = "r", chunk_size=int(10e5)
    ):
        self.h5filepath = h5filepath
        self.mode = mode
        self.chunk_size = chunk_size
        self.h5_fp: h5py.File = None
        self.chrom_container_cache = {}
        self.log = logging.getLogger("NET:MetH5")
        self.h5_fp = h5py.File(self.h5filepath, mode=self.mode)

    def __enter__(self):
        return self

    def close(self):
        self.h5_fp.close()

    def __exit__(self, exittype, exitvalue, traceback):
        self.close()

    def create_or_extend(self, parent_group, name, shape, data, **kwargs):
        if name not in parent_group.keys():
            parent_group.create_dataset(name=name, shape=shape, data=data, **kwargs)
        else:
            num_data = data.shape[0] if hasattr(data, "shape") else len(data)
            ds = parent_group[name]
            old_shape = ds.shape
            new_shape = (
                old_shape[i] + (num_data if i == 0 else 0)
                for i in range(len(old_shape))
            )
            ds.resize(new_shape)

            ds[old_shape[0] :] = data

            self.log.debug("Extended from %s to %s" % (old_shape, ds.shape))

    def add_to_h5_file(self, cur_df):
        main_group = self.h5_fp.require_group("chromosomes")

        for chrom in set(cur_df["chromosome"]):
            self.log.debug("Adding sites from chromosome %s to h5 file" % chrom)

            chrom_calls = cur_df.loc[cur_df["chromosome"] == chrom]
            n = chrom_calls.shape[0]
            read_names = [read.encode() for read in chrom_calls["read_name"]]
            read_name_len = len(read_names[0])
            assert all([len(read) for read in read_names])

            chrom_group = main_group.require_group(chrom)

            self.create_or_extend(
                parent_group=chrom_group,
                name="range",
                shape=(n, 2),
                dtype=int,
                data=chrom_calls[["start", "end"]],
                compression="gzip",
                chunks=(self.chunk_size, 2),
                maxshape=(None, 2),
            )
            self.create_or_extend(
                parent_group=chrom_group,
                name="llr",
                shape=(n,),
                dtype=float,
                data=chrom_calls["log_lik_ratio"],
                compression="gzip",
                chunks=(self.chunk_size,),
                maxshape=(None,),
            )
            self.create_or_extend(
                parent_group=chrom_group,
                name="read_name",
                shape=(n,),
                dtype="S%d" % read_name_len,
                data=read_names,
                compression="gzip",
                chunks=(self.chunk_size,),
                maxshape=(None,),
            )

            # TODO think of a way to do this that doesn't require loading everything
            # in memory
            sort_order = argsort_ranges(chrom_group)
            logging.debug("Re-sorting h5 entries for chromosome %s" % chrom)
            chrom_group["range"][:] = np.array(chrom_group["range"])[sort_order]
            chrom_group["llr"][:] = np.array(chrom_group["llr"])[sort_order]
            chrom_group["read_name"][:] = np.array(chrom_group["read_name"])[sort_order]

    def parse_and_add_nanopolish_file(self, nanopolish_file: Union[str, Path]):
        cur_df = pd.read_csv(nanopolish_file, sep="\t", dtype={"chromosome": str})
        self.add_to_h5_file(cur_df)

    def get_chromosomes(self) -> List[str]:
        return [str(k) for k in self.h5_fp["chromosomes"].keys()]

    def __getitem__(self, chromosome: str):
        if chromosome not in self.h5_fp["chromosomes"].keys():
            return None
        if chromosome in self.chrom_container_cache.keys():
            return self.chrom_container_cache[chromosome]
        else:
            ret = ChromosomeContainer(
                self.h5_fp["chromosomes"][chromosome], self.chunk_size
            )
            self.chrom_container_cache[chromosome] = ret
            return ret

    def create_chunk_index(self, force_update=False):
        for chromosome in self.get_chromosomes():
            self[chromosome].create_chunk_index(force_update=force_update)

    def annotate_read_groups(
        self, read_group_key: str, map: Dict[str, int], exists_ok=False, overwrite=False
    ):
        for chromosome in self.get_chromosomes():
            chr_g = self.h5_fp["chromosomes"][chromosome]
            rg_g = chr_g.require_group("read_groups")
            if read_group_key in rg_g.keys():
                if not exists_ok:
                    raise ValueError(
                        "Cannot annotate read groups - group assignment with this key "
                        "already exists"
                    )
                elif not overwrite:
                    continue

            rg_assignment = [
                map.get(read.decode(), -1) for read in chr_g["read_name"][:]
            ]
            rg_ds = rg_g.require_dataset(
                name=read_group_key,
                dtype=int,
                shape=(len(rg_assignment),),
                maxshape=(None,),
            )
            rg_ds[:] = rg_assignment
