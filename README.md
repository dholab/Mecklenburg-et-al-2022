# Novel Virus Discovery Pipeline & Consensus Sequence for Mecklenburg et al. 2022

This repository contains two files: 1) a snakefile script, which was used with Snakemake workflow manager, and 2) an annotated consensus sequence for the novel Hepatitis A virus variant documented in Mecklenberg et al. 2022, _A New Variant of Hepatitis A Virus Causing Transient Liver Enzyme Elevations in Mauritius-origin Laboratory-housed Cynomolhus Macaques_.

## Novel Virus Discovery Pipeline

The novel virus discovery (NVD) pipeline was run with compute resources from the Center for High Throughput Computing at the University of Wisconsin-Madison. These resources included a Docker container with all the pipeline's software dependencies, which is available at the registry `dockerreg.chtc.wisc.edu/dhoconno/nvd:26164b`.

Briefly, the pipeline goes through a number of steps to determine whether a previously undescribed virus is present in the input sequencing reads:

1. The pipeline uses many scripts from the [BBTools suite](http://sourceforge.net/projects/bbmap/), with the first being bbmask.sh. This script is used to mask host-specific reads
2. Next, the pipeline repairs paired-ed reads with repair.sh
3. The pipeline then trims Illumina adapters from the sequencing process with bbduk.sh.
4. PhiX, human, and broad array of metagenomic reads are then removed with bbmap.sh
5. Reads are then deduplicated with dedupe.sh
6. Paired and unpaired reads are then merged with bbmerge.sh.
7. Merged, trimmed reads are then assembled de novo with SPAdes (DOI: [https://doi.org/10.1002/cpbi.102](https://doi.org/10.1002/cpbi.102)).
8. Short, low-complexity contigs are removed from the SPAdes assemblies with bbmask.sh and reformat.sh
9. The pipeline then classifies contigs with megablast.
10. Unclassifiable contigs representing novel viruses are then classified with blastn.

### Running the workflow

To run the pipeline, we submitted a job for each sample with an HTCondor submit file. Each job was supplied with 8 CPUs, 64 gigabytes of RAM, and 250 gigabytes of disk space and was run in the Docker container `dockerreg.chtc.wisc.edu/dhoconno/nvd:26164b`. If you are using this pipeline outside the context of the HTCondor workfload manager, we recommend you launch the container manually with:

```
docker run \
--user $(id -u):$(id -g) -it -v $(pwd):/scratch -w /scratch \
dockerreg.chtc.wisc.edu/dhoconno/nvd:26164b /bin/bash
```

Next, we specify parameters for the pipeline. Like with compute resources for the pipeline, we specified parameters in an HTCondor submit file, though they can also be specified as environment variables in the command line. Below is how we specified the parameters for the discovery of this virus:

```
# READ_TYPE must equal se (single-end FASTQ) or pe (paired-end FASTQ)
READ_TYPE = pe

# EXPERIMENT is OC lab experiment number
# if not in OC lab, set to arbitrary identifer
EXPERIMENT = 26401

# HOST_GENOME_DB is path to reference genome
HOST_GENOME_DB = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz

# BLAST_NT is path to NCBI BLAST NT database tarball
BLAST_NT = /staging/groups/oconnor_group/nvd/blast_20210603.tar.gz

# DEPLETION_DB is path to tarball with depletion databases
DEPLETION_DB = /staging/groups/oconnor_group/nvd/21187_depletion.tar.gz

# READ_PATH is path to reads on execute node
READ_PATH = .

# OUT_DIR is path to folder on staging server where output should be saved
# CHTC wants large files on staging server, not submit server
OUT_DIR = /staging/groups/oconnor_group/26401/

# SNAKEFILE is path to snakefile script
SNAKEFILE = NVD_pipeline_snakefile.py
```

Next, we create a comma-delimited file with 3 headerless columns: 1) sample name, 2) absolute path to a read 2 fastq, and 3) absolute path to a read 1 fastq. In our HTCondor submit file, each row is then used to queue a job with the final line of the file: `queue SAMPLE_NAME,R1_FASTQ,R2_FASTQ from fastq_list.txt`. If you are running this without HTCondor, we recommend you set the file paths and sample names as environment variables, like the parameters above.

With the docker container launched, parameters set, and paths to the input files configured, you are ready to invoke the workflow with the following command:

```
snakemake --snakefile $(SNAKEFILE) --verbose -p --cores $(CPUS) \
--config read_type=$(READ_TYPE) \
experiment=$(EXPERIMENT) \
depletion_db=$(DEPLETION_DB) \
blast_nt_db=$(BLAST_NT) \
read_folder=$(READ_PATH) \
host_specific_ref_genome=$(HOST_GENOME_DB) \
r1_fastq=$(R1_FASTQ) \
r2_fastq=$(R2_FASTQ) \
out_dir=$(OUT_DIR) \
sample_name=$(SAMPLE_NAME)
```

## Consensus Sequence

After the initial results of the NVD pipeline indicated the presence of a novel, Hepatitis A Virus variant, we reprocessed reads manually in [Geneious Prime 2021.2.1](https://www.geneious.com). Like in the NVD pipeline, we first trimmed adapters, low quality sequences, and 19bp from 5’ end of each sequence. We then produced synthetic long reads by merging the trimmed reads with bbmerge, and then mapped these reads to NCBI RefSeq # NC_001489, a Hepatitis A virus reference sequence, with the Geneious mapper. To ensure we detected novel virus reads, we required at least 100 bases of overlap between the reads and the reference sequence and tolerated up to 40% base mismatches. We then generated a prototype consensus sequence from the sample with the best coverage, corrected it by mapping reads from that sample onto its consensus using Multiple Sequence Comparison by Log-Expectation align (MUSCLE; Edgar, 2004), transferred annotations from the human HAV NC_001489, trimmed low quality ends, and translated the consensus sequence into amino acids. Finally, we ran protein and nucleotide BLAST on the 7490 bp consensus sequence of the newly identified HAV, which resulted in several matches that were all approximately 20% nucleotide divergent.

The final, prototypical consensus sequence is available here as `MueHAV-consensus.gb`.
