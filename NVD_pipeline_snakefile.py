import glob
import os
import shutil
import tarfile
import http.client
import urllib
import time
import pandas as pd
from Bio import SeqIO
import sys
import re
import more_itertools
import subprocess
import labkey
from labkey.api_wrapper import APIWrapper
from labkey.query import Pagination, QueryFilter

# Arguments from CHTC --------------------------------------------------------------------------------

READ_TYPE = str(config['read_type'])
EXPERIMENT = str(config['experiment'])
BBMAP_DEPLETION_SOURCE = config['depletion_db']
BLAST_NT_SOURCE = config['blast_nt_db']
READS_FOLDER = config['read_folder']
HOST_SPECIFIC_REF_GENOME = config['host_specific_ref_genome']
R1_FASTQ = config['r1_fastq']
R2_FASTQ = config['r2_fastq']
OUTPUT_DIR = config['out_dir']
SAMPLE_NAME = config['sample_name']


# get name of BLAST folder
def rchop(s, sub):
	return s[:-len(sub)] if s.endswith(sub) else s

# names of files after copying to execute node
BLAST_DB_NAME = os.path.basename(rchop(BLAST_NT_SOURCE, '.tar.gz'))

R1_FASTQ_EXECUTE_NODE = 'input_S1_L001_R1_001.fastq.gz'
R2_FASTQ_EXECUTE_NODE = 'input_S1_L001_R2_001.fastq.gz'

# test sanity of read type
# only allowed values are pe (paired-end) or se (single-end)
if READ_TYPE != 'se' and READ_TYPE != 'pe':
	print('Read type needs to be specified as "pe" or "se". Aborting now.')
	sys.exit()

# path to read cleaning databases
READ_CLEAN_DB_PATH = 'ref'

# Rules ----------------------------------------------------------------------------------

onerror:
	# handler for error
	# mainly so CHTC returns .tar.gz that can be inspected to understand why snakemake failed

	# make compressed tarball when processing completes
	with tarfile.open(SAMPLE_NAME + '-error.tar.gz', "w:gz") as tar:
		tar.add(SAMPLE_NAME + '/', arcname=os.path.basename(SAMPLE_NAME + '/'))
		
	# copy output to staging server
	shell('cp ' + SAMPLE_NAME + '-error.tar.gz ' + OUTPUT_DIR)
	
	# remove compressed tarball so it doesn't also transfer to submit server
	shell('rm ' + SAMPLE_NAME + '-error.tar.gz ')
	

rule all:
	'''require processing through LabKey upload and mapping of cleaned BAM reads against SPAdes contigs'''
	# snakemake uses rule all to specify destination files
	# the wildcards expand the sample names derived from FASTQ files in the READS_FOLDER into the necessary output files
	input:
		SAMPLE_NAME + '/blast/abridged.classified_contigs.txt',
		SAMPLE_NAME + '/labkey/' + SAMPLE_NAME + '.labkey.tkn',
		SAMPLE_NAME + '/blast/hq.virus.contigs.txt',
		SAMPLE_NAME + '/virus_bam/' + SAMPLE_NAME + '.virus_sorted.bam'
	run:
		# remove empty directories
		shell('find ' + ' -type d -empty -delete')

		# make compressed tarball when processing completes
		with tarfile.open(SAMPLE_NAME + '.tar.gz', "w:gz") as tar:
			tar.add(SAMPLE_NAME  + '/', arcname=os.path.basename(SAMPLE_NAME + '/'))

		# copy results to CHTC staging server
		# copy to output_destination
		shell('mkdir -p ' + OUTPUT_DIR + ' && cp ' + SAMPLE_NAME  + '.tar.gz ' + OUTPUT_DIR)
		
		# remove compressed tarball so it doesn't also transfer to submit server
		shell('rm ' + SAMPLE_NAME  + '.tar.gz ')
		
		# then remove results folder
		shutil.rmtree(SAMPLE_NAME  + '/')

		# remove ref databases from working directory
		# important because BLAST databases are 120GB
		# disabled because it doesn't transfer back to CHTC submit node and it increases time for testing on dhogal2
		# shutil.rmtree('ref/')
		# remove copied depletion and BLAST databases so they don't transfer back to execute node
		shell('rm depletion.tar.gz nt.tar.gz')

		# print arguments to snakemake log
		# important to know exactly what versions of databases were used for analysis
		print('--Arguments--')
		print('Experiment number: ' + EXPERIMENT)
		print('bbmap depletion and trim databases: ' + BBMAP_DEPLETION_SOURCE)
		print('BLAST nt database: ' + BLAST_NT_SOURCE)

## Copy files from CHTC transfer server to execute node ##

rule copy_fastq:
	'''
	copy FASTQ files from CHTC transfer server to execute node
	'''
	input:
		R1_FASTQ,
		R2_FASTQ
	output:
		temp(R1_FASTQ_EXECUTE_NODE),
		temp(R2_FASTQ_EXECUTE_NODE)
	threads: 1
	run:
		shell('cp {input[0]} {output[0]}')
		shell('cp {input[1]} {output[1]}')

rule decompress_depletion_db:
	'''
	copy the depletion .tar.gz archive to the working directory
	decompress the depletion .tar.gz archive in the working directory
	'''
	output:
		temp('depletion_decompressed.txt')
	threads: 1
	run:
		# if depletion database has not been decompressed (based on not having a reference director), do so here
		shell('if [ ! -f "ref/adapters.fa" ] ; then cp ' + BBMAP_DEPLETION_SOURCE + ' depletion.tar.gz && tar -xzf depletion.tar.gz && touch {output[0]} ; fi')
	
rule decompress_blast_db:
	'''
	decompress BLAST .tar.gz archive once per execute node
	create temporary output file to indicate that decompressed BLAST folder is available
	'''
	input:
		BLAST_NT_SOURCE
	output:
		temp('blast_decompressed.txt')
	threads: workflow.cores
	run:
		# copy BLAST database .tar.gz from staging to execute node
		shell('cp ' + BLAST_NT_SOURCE + ' nt.tar.gz')
		
		# use pigz if it is available
		shell('if [ -x "$(command -v pigz)" ] ; \
			then pigz -dc nt.tar.gz | tar -xf -  && touch {output[0]}; \
			else tar -xzf nt.tar.gz  && touch {output[0]} ; fi ;')	

## Prepare sequencing reads for de novo assembly ##

rule mask_genome:
	'''create low-complexity masked reference genome to remove host-specific reads
	necessary since human read removal isn't effective for divergent NHP'''
	params:
		FN = HOST_SPECIFIC_REF_GENOME
	output:
		temp('ref/host-genome.fa.gz'),
		temp('ref/host-genome.masked.fa.gz')
	threads: workflow.cores
	run:
		shell('if [ ! -f {output[0]} ] ; then wget -O {output[0]} {params.FN} ; fi')
		shell('bbmask.sh -Xmx64g in={output[0]} out={output[1]}')

if READ_TYPE == 'pe':
	rule repair_fastq:
		'''
		run bbamp repair.sh on input files to make sure that they are properly paired
		'''
		input:
			R1_FASTQ_EXECUTE_NODE,
			R2_FASTQ_EXECUTE_NODE
		output:
			temp(SAMPLE_NAME + '/repaired/' + SAMPLE_NAME + '.fastq.gz')
		threads: 4
		run:
			shell('repair.sh \
			in={input[0]} \
			in2={input[1]}\
			out={output[0]}')
		
	rule trim_adapters:
		"""run bbduk to remove adapters. use adapter definitions from bbmap.
		use parameters from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/"""

		input:
			SAMPLE_NAME + '/repaired/' + SAMPLE_NAME + '.fastq.gz',
			'depletion_decompressed.txt'
		output:
			temp(SAMPLE_NAME + '/adapter_trimming/' + SAMPLE_NAME + '.fastq.gz')
		log:
			SAMPLE_NAME + '/logs/adapter_trimming/' + SAMPLE_NAME + '.log'
		threads: workflow.cores
		run:
			# remove adapters
			shell('bbduk.sh \
			-Xmx64g \
			in={input[0]} \
			out={output} \
			ref=' + READ_CLEAN_DB_PATH + '/adapters.fa \
			threads={threads} \
			ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> {log}')

if READ_TYPE == 'se':
	rule trim_adapters:
		"""run bbduk to remove adapters. use adapter definitions from bbmap.
		use parameters from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/"""

		input:
			R1_FASTQ_EXECUTE_NODE,
			'depletion_decompressed.txt'
		output:
			temp(SAMPLE_NAME + '/adapter_trimming/' + SAMPLE_NAME + '.fastq.gz')
		log:
			SAMPLE_NAME + '/logs/adapter_trimming/' + SAMPLE_NAME + '.log'
		threads: workflow.cores
		run:
			# remove adapters
			shell('bbduk.sh \
						-Xmx64g \
						in={input[0]} \
						out={output} \
						ref=' + READ_CLEAN_DB_PATH + '/adapters.fa \
						threads={threads} \
						ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> {log}')

rule remove_host_specific:
	"""run host-specific reads"""
	input:
		SAMPLE_NAME + '/adapter_trimming/' + SAMPLE_NAME + '.fastq.gz',
		'ref/host-genome.masked.fa.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_host_specific.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '.remove_host_specific.log'
	threads: workflow.cores
	run:
		# remove host_specific reads
		shell('bbmap.sh \
		-Xmx64g \
		usemodulo=t \
		int=t \
		in={input[0]} \
		outu={output} \
		ref={input[1]} \
		threads={threads} \
		minid=0.85 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 trimq=10 untrim ow=t 2> {log}')

rule remove_phix:
	"""run bbduk to remove phiX reads."""
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_host_specific.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_phix.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '.remove_phix.log'
	threads: workflow.cores
	run:
		# remove phix
		shell('bbduk.sh \
		-Xmx64g \
		in={input} \
		outu={output} \
		path=' + READ_CLEAN_DB_PATH + '/phix_adapters.fa.gz \
		threads={threads} \
		k=31 hdist=1 ow=t 2> {log}')

rule remove_human:
	"""remove human sequences
	   make human database with:
	   bbmap.sh ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz path=./hg19_main_mask_ribo_animal_allplant_allfungus -Xmx64g -usemodulo=t
	   usemodulo to reduce RAM usage and enable 16Gb CHTC memory allocations"""
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_phix.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_human.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '.remove_human.log'
	threads: workflow.cores
	run:
		shell('bbmap.sh \
		-Xmx64g \
		usemodulo=t \
		int=t \
		in={input} \
		outu={output} \
		path=' + READ_CLEAN_DB_PATH + '/hg19_main_mask_ribo_animal_allplant_allfungus \
		threads={threads} \
		minid=0.85 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 trimq=10 untrim ow=t 2> {log}')

rule remove_broad_metagenomics:
	"""remove Broad metagenomics contaminants
	   from https://platform.dnanexus.com/projects/F8PQ6380xf5bK0Qk0YPjB17P/data/resources/production_defaults"""
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_human.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_broad_metagenomics.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '.remove_broad_metagenomics.log'
	threads: workflow.cores
	run:
		shell('bbmap.sh \
		-Xmx64g \
		int=t \
		in={input} \
		outu={output} \
		path=' + READ_CLEAN_DB_PATH + '/metagenomics_contaminants_v3 \
		threads={threads} \
		minid=0.85 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 trimq=10 untrim ow=t 2> {log}')

rule remove_fusedEPmasked2:
	"""remove bbmap fusedEPmasked2 bacterial contaminants
	   from http://seqanswers.com/forums/archive/index.php/t-42552.html"""
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_broad_metagenomics.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_fusedEPmasked2.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '.remove_fusedEPmasked2.log'
	threads: workflow.cores
	run:
		shell('bbmap.sh \
		-Xmx64g \
		int=t \
		in={input} \
		outu={output} \
		path=' + READ_CLEAN_DB_PATH + '/fusedEPmasked2 \
		threads={threads} \
		minid=0.85 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 trimq=10 untrim ow=t 2> {log}')

rule remove_fusedERPBBmasked2:
	"""remove bbmap fusedEPmasked2 bacterial contaminants
	   from http://seqanswers.com/forums/archive/index.php/t-42552.html"""
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_fusedEPmasked2.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_fusedERPBBmasked2.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '.remove_fusedERPBBmasked2.log'
	threads: workflow.cores
	run:
		shell('bbmap.sh \
		-Xmx64g \
		int=t \
		in={input} \
		outu={output} \
		path=' + READ_CLEAN_DB_PATH + '/fusedERPBBmasked2 \
		threads={threads} \
		minid=0.85 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 trimq=10 untrim ow=t 2> {log}')

rule remove_metag_v3_ncRNA_mRNA_mitRNA_consensus:
	"""remove Broad metag_v3_ncRNA_mRNA_mitRNA_consensus
	   from https://platform.dnanexus.com/projects/F8PQ6380xf5bK0Qk0YPjB17P/data/resources/production_defaults"""
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_fusedERPBBmasked2.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_metag_v3_ncRNA_mRNA_mitRNA_consensus.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/remove_contaminants/' + SAMPLE_NAME + '_remove_metag_v3_ncRNA_mRNA_mitRNA_consensus.log'
	threads: workflow.cores
	run:
		shell('bbmap.sh \
		-Xmx64g \
		int=t \
		in={input} \
		outu={output} \
		path=' + READ_CLEAN_DB_PATH + '/metag_v3_ncRNA_mRNA_mitRNA_consensus \
		threads={threads} \
		minid=0.85 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 trimq=10 untrim ow=t 2> {log}')

rule deduplicate_reads:
	'''run bbmap dedupe to remove duplicate sequences'''
	input:
		SAMPLE_NAME + '/remove_contaminants/' + SAMPLE_NAME + '.remove_metag_v3_ncRNA_mRNA_mitRNA_consensus.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/dedupe/' + SAMPLE_NAME + '.dedupe.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/dedupe/' + SAMPLE_NAME + '.log'
	threads: workflow.cores
	run:
		shell('dedupe.sh \
		-Xmx64g \
		int=t \
		in={input} \
		out={output} \
		threads={threads} 2> {log}')

rule merge_reads:
	'''run bbmerge and save both merged and unpaired reads'''
	input:
		SAMPLE_NAME + '/dedupe/' + SAMPLE_NAME + '.dedupe.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/merged/' + SAMPLE_NAME + '.merged.fastq.gz'),
		temp(SAMPLE_NAME + '/merged/' + SAMPLE_NAME + '.unmerged.fastq.gz')
	log:
		SAMPLE_NAME + '/logs/merged/' + SAMPLE_NAME + '.log'
	threads: workflow.cores
	run:
		shell('bbmerge.sh \
		-Xmx64g \
		int=t \
		in={input} \
		out={output[0]} \
		outu={output[1]} \
		threads={threads} 2> {log}')

## De novo assemble pre-processed reads ##

rule spades_assemble:
	'''
	run SPAdes in RNA mode
	this empirically gives the best viral RNA contigs
	though sometimes viruses are represented in 2 or more nearly identical contigs
	'''
	input:
		SAMPLE_NAME + '/merged/' + SAMPLE_NAME + '.merged.fastq.gz',
		SAMPLE_NAME + '/merged/' + SAMPLE_NAME + '.unmerged.fastq.gz'
	output:
		SAMPLE_NAME + '/spades_fasta/transcripts.fasta'
	log:
		SAMPLE_NAME + '/logs/spades_fasta/' + SAMPLE_NAME + '.log'
	threads: workflow.cores
	run:
		# de novo assemble
		# if single-end reads
		if READ_TYPE == 'se':
			shell('spades.py --rna -t {threads} -s {input[1]}  -o ' + SAMPLE_NAME + '/spades_fasta/')

		# if paired-end reads
		if READ_TYPE == 'pe':
			shell('spades.py --rna -t {threads} --12 {input[1]} -s {input[0]}  -o ' + SAMPLE_NAME + '/spades_fasta/')

		# remove intermediate files
		shell('rm -rf ' + SAMPLE_NAME + '/spades_fasta/corrected')
		shell('rm -rf ' + SAMPLE_NAME + '/spades_fasta/K127')
		shell('rm -rf ' + SAMPLE_NAME + '/spades_fasta//misc')
		shell('rm -rf ' + SAMPLE_NAME + '/spades_fasta//split_input')
		shell('rm -rf ' + SAMPLE_NAME + '/spades_fasta//tmp')
		shell('rm ' + SAMPLE_NAME + '/spades_fasta/before_rr.fasta')
		shell('rm ' + SAMPLE_NAME + '/spades_fasta/dataset.info')
		shell('rm ' + SAMPLE_NAME + '/spades_fasta/hard_filtered_transcripts.fasta')
		shell('rm ' + SAMPLE_NAME + '/spades_fasta/input_dataset.yaml')
		shell('rm ' + SAMPLE_NAME + '/spades_fasta/soft_filtered_transcripts.fasta')

		# copy log file generated by SPAdes to snakemake log
		shell('cp ' + SAMPLE_NAME + '/spades_fasta/spades.log {log}')

rule remove_short_contigs:
	'''remove contigs shorter than 500bp as these are likely not going to informative. They
	also clog the blast classification pipeline.'''
	input:
		SAMPLE_NAME + '/spades_fasta/transcripts.fasta'
	output:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp.fasta'
	log:
		SAMPLE_NAME + '/logs/spades_fasta/' + SAMPLE_NAME + '_500bp.log'
	threads: workflow.cores
	run:
		shell('reformat.sh \
		-Xmx64g \
		minlength=500 \
		in={input} \
		out={output} \
		2> {log}')

rule mask_low_complexity_short_contigs:
	'''use bbmask to N-mask low complexity sequence to reduce BLAST runtimes.
	remove N from ends of sequences by piping to reformat.sh
	require at least 100bp of called bases after trimming
	necessary becuase otherwise BLAST can hang'''
	input:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp.fasta'
	output:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp_masked.fasta'
	log:
		SAMPLE_NAME + '/logs/spades_fasta/' + SAMPLE_NAME + '_500bp_masked.log'
	threads: workflow.cores
	run:
		shell('bbmask.sh \
		-Xmx64g \
		in={input} \
		out=stdout.fa \
		| reformat.sh \
		in=stdin.fa \
		out={output} \
		qtrim=t \
		minconsecutivebases=100 \
		2> {log}')

## Classify reads with BLAST ##

rule megablast_classify:
	'''
	run megablast on spades assembly output
	pass unclassified sequences to blastn
	'''
	input:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp_masked.fasta',
		'blast_decompressed.txt'
	output:
		SAMPLE_NAME + '/blast/megablast.out.txt'
	threads: workflow.cores
	run:
		shell('blastn -task megablast \
		-db ' + BLAST_DB_NAME + '/nt \
		-query {input[0]} \
		-num_threads {threads} \
		-outfmt "6 qseqid qlen sseqid stitle length pident evalue bitscore" \
		-out {output[0]}')
		
rule annotate_megablast_results:
	'''
	create file with megablast results annotated with megablast task
	'''
	input:
		SAMPLE_NAME + '/blast/megablast.out.txt'
	output:
		SAMPLE_NAME + '/blast/megablast.annotated.out.txt'
	threads: 1
	run:
		# iterate over BLAST results, add task and sample name, and write to merged file
		with open(input[0], 'r') as f:
			lines = f.readlines()
		lines = ['megablast\t' + SAMPLE_NAME + '/\t' +line for line in lines]
		with open(output[0], 'w') as f:
			f.writelines(lines)

rule remove_megablast_mapped_contigs:
	'''
	create list of contigs that are not classified by megablast
	used as input for more sensitive blastn searching
	'''
	input:
		SAMPLE_NAME + '/blast/megablast.out.txt',
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp_masked.fasta',
	output:
		SAMPLE_NAME + '/blast/megablast.classified.txt',
		SAMPLE_NAME + '/blast/megablast.pruned.fa'
	threads: 1
	run:
		# save unique IDs in BLAST results table to a python set
		# get IDs from first column of BLAST output table

		unique_ids = set()
		with open(input[0]) as f:
			for line in f:
				unique_ids.add(line.rstrip('\n').split('\t')[0] + '\n')

		# write unique_ids to file
		with open(output[0],'w') as f:
			for i in unique_ids:
				f.write(i)

		# use bbmap filterbyname.sh to generate output FASTA with unmapped contigs
		shell('filterbyname.sh \
		in={input[1]} \
		out={output[1]} \
		names={output[0]}')
		
rule blastn_classify:
	'''
	run blastn on sequences not classified by megablast
	'''
	input:
		SAMPLE_NAME + '/blast/megablast.pruned.fa',
		'blast_decompressed.txt'
	output:
		SAMPLE_NAME + '/blast/blastn.out.txt'
	threads: workflow.cores
	run:
		shell('blastn -task blastn \
		-db ' + BLAST_DB_NAME + '/nt \
		-query {input[0]} \
		-num_threads {threads} \
		-outfmt "6 qseqid qlen sseqid stitle length pident evalue bitscore" \
		-out {output[0]}')

rule annotate_blastn_results:
	'''
	create file with blastn results annotated with blastn task
	'''
	input:
		SAMPLE_NAME + '/blast/blastn.out.txt'
	output:
		SAMPLE_NAME + '/blast/blastn.annotated.out.txt'
	threads: 1
	run:
		# iterate over BLAST results, add task and sample name, and write to merged file
		with open(input[0], 'r') as f:
			lines = f.readlines()
		lines = ['blastn\t' + SAMPLE_NAME + '/\t' +line for line in lines]
		with open(output[0], 'w') as f:
			f.writelines(lines)

rule merge_annotated_blast_results:
	'''
	create file with merged megablast and blastn results
	'''
	input:
		SAMPLE_NAME + '/blast/megablast.annotated.out.txt',
		SAMPLE_NAME + '/blast/blastn.annotated.out.txt'
	output:
		SAMPLE_NAME + '/blast/merged.out.txt'
	threads: 1
	run:
		shell('cat {input[0]} {input[1]} > {output[0]}')
		
rule abridge_blast_results:
	'''
	create file with only top 5 blast hits per query sequence by bit score
	'''
	input:
		SAMPLE_NAME + '/blast/merged.out.txt'
	output:
		SAMPLE_NAME + '/blast/abridged.classified_contigs.txt'
	run:
		#expects blast_results will be in the format of merged_blast_output with two leading columns showing blast_task and sample_id
		if os.path.getsize(input[0]) > 0:
			df = pd.read_csv(input[0], sep='\t',header=None, quoting=3)
			df.sort_values([2, 9], ascending=[True, False]).groupby(2).head(5).reset_index(drop=True).to_csv(output[0], sep='\t')

rule post_process_blast:
	'''evaluate BLAST results and identify contigs that are confidently classified as viral hits
	make a CSV files of these BLAST hits and save sequences as FASTA files.
	also save FASTA files for spades contigs that are not classified by BLAST and may be truly novel'''
	input:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp_masked.fasta',
		SAMPLE_NAME + '/blast/abridged.classified_contigs.txt'
	output:
		SAMPLE_NAME + '/blast/hq.virus.contigs.txt',
		SAMPLE_NAME + '/blast/hq.virus.contigs.fasta',
		SAMPLE_NAME + '/blast/unclassified.contigs.fasta'
	log:
		SAMPLE_NAME + '/logs/blast/postprocess.log'
	threads: 1
	run:
		# create output files in case input file is empty or there are no virus hits in blast results
		shell('touch {output[0]} && touch {output[1]}')

		def create_blast_dataframe(blast_results):
			'''import CSV of BLAST results to pandas dataframe.
			return dataframe with proper column names'''

			# import abridged BLAST classification
			BLAST_CLASSIFICATIONS_DF = pd.read_csv(blast_results, sep='\t', index_col=0)

			# rename columns
			BLAST_CLASSIFICATIONS_DF.rename(columns={'0': 'algorithm',
			'1': 'sample',
			'2': 'qseqid',
			'3': 'qlen',
			'4': 'sseqid',
			'5': 'stitle',
			'6': 'length',
			'7': 'pident',
			'8': 'evalue',
			'9': 'bitscore'
			}, inplace=True)

			return BLAST_CLASSIFICATIONS_DF

		def get_hq_virus_contigs(blast_df, fasta_contigs):
			'''identify contigs where all BLAST hits are viruses
			input is blast result dataframe'''

			# filter on column 5 descrption
			# save entries with term virus
			VIRUS_CLASSIFICATIONS_DF = blast_df[blast_df['stitle'].str.contains('virus')]

			# if dataframe is empty (meaning no virus hits)
			# quit script
			if VIRUS_CLASSIFICATIONS_DF.empty:
				print('No virus contigs in BLAST results')
				return

			# count number of BLAST entries for each contig
			# though most should have 5, if a contig is a highly divergent virus there might be less than 5
			BLAST_HITS_PER_NODE_DF = blast_df.groupby(['qseqid']).size()

			# count number of viral entries for each contig
			VIRUS_HITS_PER_NODE_DF = VIRUS_CLASSIFICATIONS_DF.groupby(['qseqid']).size()

			# make dataframe with total BLAST hits in 0 and viral BLAST hits in 1
			COMPARE_VIRUS_CONTIGS_DF = pd.concat([BLAST_HITS_PER_NODE_DF, VIRUS_HITS_PER_NODE_DF], axis=1).reset_index()

			# rename columns
			COMPARE_VIRUS_CONTIGS_DF.rename(columns={'index': 'qseqid',
			0: 'blast_ct',
			1: 'virus_ct'
			}, inplace=True)

			# convert NA values to 0 in column '1'
			COMPARE_VIRUS_CONTIGS_DF.fillna(0)

			# filter on columns where values are equal in 0 and 1
			# these are columns where all available BLAST hits in abrdiged report are from a virus
			HQ_VIRUS_CONTIGS = COMPARE_VIRUS_CONTIGS_DF[COMPARE_VIRUS_CONTIGS_DF['blast_ct'] == COMPARE_VIRUS_CONTIGS_DF['virus_ct']]

			# get names of HQ virus contig nodes as list
			# column with node names is named '2' (yes, as a string)
			HQ_VIRUS_NODES = HQ_VIRUS_CONTIGS['qseqid'].tolist()

			# get rows in BLAST_CLASSIFICATIONS_DF that have node name in list HQ_VIRUS_NODES
			# these are the high confidence virus sequences that should be examined further
			# export as CSV
			HQ_VIRUS_TABLE = blast_df[blast_df['qseqid'].isin(HQ_VIRUS_NODES)]
			HQ_VIRUS_TABLE.to_csv(output[0], sep='\t')

			# save FASTA files for nodes that are HQ virus contigs
			# use bbmap filterbyname.sh
			HQ_VIRUS_NODES_NAMES = ','.join(HQ_VIRUS_NODES)

			# default behavior is to exclude matching reads
			# we want to specifically include matching sequences
			cmd = ['filterbyname.sh',
				  'include=t',
				  'in=' + fasta_contigs,
				  'names=' + HQ_VIRUS_NODES_NAMES,
				  'out=' + SAMPLE_NAME + '/blast/hq.virus.contigs.fasta',
				  'ow=t']

			subprocess.call(cmd)

			return None

		def get_unclassified_contigs(blast_df, fasta_contigs):
			'''save FASTA files from interesting contigs where there are NO BLAST matches whatsoever'''

			# get all sequence names from spades contigs
			SPADES_CONTIGS = []
			for record in SeqIO.parse(fasta_contigs, "fasta"):
				SPADES_CONTIGS.append(record.id)

			# convert sequence names to pandas series
			SPADES_NODES = pd.Series(v for v in SPADES_CONTIGS)

			# get unique nodes from BLAST hits
			BLAST_NODES = blast_df['qseqid'].unique()

			# get sequences whose names are in SPADES_CONTIGS but not BLAST_NODES using set logic
			UNCLASSIFIED_NODES = set(SPADES_NODES) - set(BLAST_NODES)

			# if UNCLASSIFIED_NODES is not empty, extract FASTA sequences of these nodes
			if len(UNCLASSIFIED_NODES) > 0:
				UNCLASSIFIED_NODE_NAMES = ','.join(UNCLASSIFIED_NODES)

				# default behavior is to exclude matching reads
				# we want to specifically include matching sequences
				cmd = ['filterbyname.sh',
					  'include=t',
					  'in=' + fasta_contigs,
					  'names=' + UNCLASSIFIED_NODE_NAMES,
					  'out=' + output[2]]

				subprocess.call(cmd)
			else:
				shell('touch {output[2]}')

				return None

		# make BLAST dataframe
		BLAST_DF = create_blast_dataframe(input[1])

		# get HQ virus contigs
		get_hq_virus_contigs(BLAST_DF, str(input[0]))

		# get unclassified contigs
		get_unclassified_contigs(BLAST_DF, str(input[0]))

## Process BLAST results and load into LabKey ##

rule fasta_singleline:
	"""convert multiline FASTA sequences to single line FASTA sequences. Necessary for LabKey FASTA import"""
	input:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp_masked.fasta'
	output:
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp.singleline.fasta'
	threads: 1
	shell: "awk '/^>/ {{printf(\"\\n%s\\n\",$0);next; }} {{ printf(\"%s\",$0);}}  END {{printf(\"\\n\");}}' < {input} > {output}"

rule labkey_upload:
	'''run script that uploads classified contig information and FASTA spades contigs to labkey.
	authenticates to Labkey via API key that is only able to interact with tables in this folder.
	make sure labkey and more_itertools modules installed via pip before running this rule'''
	input:
		SAMPLE_NAME + '/blast/abridged.classified_contigs.txt',
		SAMPLE_NAME + '/spades_fasta/transcripts_500bp_masked.fasta'
	output:
		SAMPLE_NAME + '/labkey/' + SAMPLE_NAME + '.labkey.tkn' # empty placeholder token to invoke this rule
	threads: 1
	log:
		SAMPLE_NAME + '/logs/labkey/' + SAMPLE_NAME + '.log'
	run:
		shell('touch {output}')

		# functions for use in conjunction with de novo assembly and BLAST classification
		# run this in a subprocess

		def insertBlastResults(experiment, blast_results):
			"""insert records into Labkey https://dholk.primate.wisc.edu/list/dho/projects/nvd/novel_virus_sequencing/"""

			# parse blast_results into python dictionary
			if os.path.getsize(blast_results) > 0:
				df = pd.read_csv(blast_results, sep='\t',header=None, quoting=3)
				blast_rows = []
				for r in df.itertuples():
					if 'blast' in r[2]:
						#save results to python list that will be used to populate labkey import
						blast_rows.append(
						{'experiment': str(experiment),
						'blast_task': str(r[2]),
						'sample_id': str(r[3]),
						'qseqid': str(r[4]),
						'qlen': str(r[5]),
						'sseqid': str(r[6]),
						'stitle': str(r[7]),
						'length': str(r[8]),
						'pident': str(r[9]),
						'evalue': str(r[10]),
						'bitscore': str(r[11]),
						'blast_db_version': str(BLAST_DB_NAME),
						'nvd_version': '25880'}
						)

				# Labkey has a hard time handling large numbers of inserts a time
				# break apart fasta_rows into smaller chunks of 1000 rows each

				blast_row_chunks = list(more_itertools.chunked(blast_rows, 1000))
				for counter, i in enumerate(blast_row_chunks):
					#insert rows into labkey
					labkey_server = 'dholk.primate.wisc.edu'
					project_name = 'dho/projects/nvd/novel_virus_sequencing/'  # Project folder name
					api = APIWrapper(labkey_server, project_name, api_key='apikey|4e80e6d4f9099b6d70f6fbbffc071e46', use_ssl=True)

					api.query.insert_rows(
					schema_name = 'oconnor_external',
					query_name  = 'nvd_blast_output',
					rows = i)

					# print updates as rows processed
					number_processed = counter * 1000
					print(str(number_processed) + ' BLAST output rows added')

		def insertFastaFile(experiment, fasta_contigs):
			# extract contig_id and sequence from fasta_contigs
			# open file, use re.finditer to capture all matches
			# write results to list that can be inserted into labkey
			fasta_rows = []

			with open(fasta_contigs, 'r') as f:
				data = (f.read())
				for match in re.finditer('>(NODE.*[0-9]*)\n([ACTGN]*)', data):
					fasta_rows.append(
					{'experiment': str(experiment),
					'sample_id': SAMPLE_NAME,
					'contig_id': match.group(1),
					'contig_sequence': match.group(2),
					'nvd_version': '25880'}
					)

			# Labkey has a hard time handling large numbers of inserts a time
			# break apart fasta_rows into smaller chunks

			fasta_row_chunks = list(more_itertools.chunked(fasta_rows, 1000))
			
			for counter, i in enumerate(fasta_row_chunks):
				#insert rows into labkey
				labkey_server = 'dholk.primate.wisc.edu'
				project_name = 'dho/projects/nvd/novel_virus_sequencing/'  # Project folder name
				api = APIWrapper(labkey_server, project_name, api_key='apikey|4e80e6d4f9099b6d70f6fbbffc071e46', use_ssl=True)

				api.query.insert_rows(
				schema_name ='oconnor_external',
				query_name ='nvd_fasta_contigs',
				rows = i)

				# print updates as rows processed
				number_processed = counter * 1000
				print(str(number_processed) + ' FASTA contigs added')

		# run labkey insert functions

		insertBlastResults(EXPERIMENT, str(input[0]))
		insertFastaFile(EXPERIMENT, str(input[1]))

## Map reads to virus contigs ##

rule map_to_virus:
	'''map all input reads to high quality virus FASTA contigs and save output as BAM files.
	this will make it easier to inspect authenticity of putative virus contigs'''
	input:
		SAMPLE_NAME + '/blast/hq.virus.contigs.fasta',
		SAMPLE_NAME + '/adapter_trimming/' + SAMPLE_NAME + '.fastq.gz'
	output:
		temp(SAMPLE_NAME + '/virus_bam/' + SAMPLE_NAME + '.virus.bam')
	log:
		SAMPLE_NAME + '/logs/virus_bam/' + SAMPLE_NAME + '.log'
	threads: workflow.cores
	run:
		# map input reads to reference files
		# save only mapped reads
		# save all ambiguous mapping results
		# allow indels up to 100bp - default can cause BAM files that do not parse in Geneious correctly
		# handle case where reference file does not contain any sequences - this causes bbmap to error
		# create empty output file if reference file is empty so next step won't fail
		shell('if [ -s {input[0]} ] ; \
			  then \
				bbmap.sh \
				-Xmx64g \
				in={input[1]} \
				ref={input[0]} \
				threads={threads} \
				minid=0.95 \
				maxindel=100 \
				outm={output[0]} \
				ambiguous=all 2> {log} ; \
			  else \
				touch {output}; \
			  fi')

		# copy FASTA file of high quality virus sequences to BAM folder
		shell('cp {input[0]} ' + SAMPLE_NAME + '/virus_bam/' + SAMPLE_NAME + '.hq.virus.contigs.fasta')

rule sort_index_virus_bam:
	'''create sorted, indexed BAM file following mapping to virus contigs. Necessary because 
	sometimes lots of virtually identical contigs are made by SPAdes and trying to load all
	mappings in a program like Geneious doesn't work well. Plus, it's better to have sorted,
	indexed BAM files'''
	input:
		SAMPLE_NAME + '/virus_bam/' + SAMPLE_NAME + '.virus.bam'
	output:
		SAMPLE_NAME + '/virus_bam/' + SAMPLE_NAME + '.virus_sorted.bam'
	threads: workflow.cores
	run:
			# handle cases where input bam is empty by creating empty output file
		shell('if [ -s {input} ] ; \
			  then \
				samtools sort -o {output} {input} && \
				samtools index {output}; \
			  else \
				touch {output}; \
			  fi')

		# handle empty files by creating if they do not exist
		shell('touch {output[0]}')