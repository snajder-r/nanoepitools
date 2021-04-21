from typing import List
import pysam
import bisect
import tqdm

class WhereDoReadsGo:
    def __init__(self, bam_files:List[str], chunk_size=int(1e4)):
        self.bam_files = bam_files
        
    def index(self):
        self.read_index = {}
        
        self.coordinate_index = {}
        
        with tqdm.tqdm(total=len(self.bam_files)) as pbar:
            for bam_file in self.bam_files:
   
                """
                More efficient to write into a fresh dict and join them later, since it may avoid collisions
                """
                batched_read_index = {}
                with pysam.AlignmentFile(bam_file, "r") as bam:
                    reference_lengths = {bam.get_reference_name(i):bam.get_reference_length(i) for i in range(bam.nreferences)}
                    for chrom, chrom_len in reference_lengths.items():
                    
                    
                    
                    for line in bam.fetch():
                        if line.is_secondary:
                            continue
                        # Read index
                        if line.query_name not in batched_read_index:
                            read_alignments = []
                            batched_read_index[line.query_name] = read_alignments
                        else:
                            read_alignments = batched_read_index[line.query_name]
                        bisect.insort(read_alignments, (line.query_alignment_start, line.query_alignment_end, line.reference_name, line.reference_start, line.reference_end))
                        
                        # Coordinate index
                        
                            chrom_coord_index = self.coordinate_index[line.reference_name]
                            
                        
                        
                self.read_index.update(batched_read_index)
                pbar.update(1)
        return self