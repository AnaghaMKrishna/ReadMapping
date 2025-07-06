# ReadMapping

This code repo maps raw reads to reference genome to obtain the alignments and based on the best consensus sequence obtained, infer phylogeny of the organisms. The application takes raw reads in FASTQ format or assembled reads in FASTA format, extracts the 16S sequences from FASTQs or FASTAs. If FASTQs are provided as input, the raw reads are mapped to muliple reference genomes using `minimap2` to obtain alignments in SAM file. The `process_sam.py` file processes each aligned read to obtain consensus bases at every position. The best consensus is decided based on the highest horizontal coverage of the reads for all positions in the reference. The assembled reads are then used with `isPCR` module to extract the amplicons, oriented in the same direction using Needleman-Wunsch global alignment algorithm implemented in `nw` module and finally aligned using Multiple Sequence Alignment(MSA) using `muscle`. The output is then used for phylogenetic inference to classify unknown prokaryotic organisms by placing them in a phylogenetic tree based on an analysis of their 16S rRNA gene sequence using `iqtree` followed by `gotree` for midpoint rooting.

All the outputs can be found under `readmapping/outputs` and the final phylogenetic tree is present in `rooted_ml_tree`.

General Usage:
- `./phyl_analysis.py [-a ASSEMBLY [ASSEMBLY ...]] -p PRIMERS [-r READS [READS ...]] [-sREF_SEQS]`
- `muscle -align <input file> -output <output alignment>`
- `iqtree2 -s <output alignment> --prefix tree/ml_tree --quiet`
- `gotree reroot midpoint -i <input treefile> -o <output rooted treefile>`

Example Usage:
- `./phyl_analysis.py -a data/assemblies/<assembly_files> -p data/primers/<16S primer file> -r data/reads/<PE read files> -s data/refs/<references> > readmapping/outputs/<aligned_output.fna>`
- `muscle -align readmapping/outputs/<aligned_output.fna> -output readmapping/outputs/<msa_output>`
- `iqtree2 -s readmapping/outputs/<msa_output> --prefix readmapping/outputs/ml_tree --quiet`
- `gotree reroot midpoint -i readmapping/outputs/ml_tree.treefile -o readmapping/outputs/rooted_ml_tree`
