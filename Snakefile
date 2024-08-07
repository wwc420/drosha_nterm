# Worflow to map small RNA sequencing data to AQseq SpikeIns first and Human genome Next
# First written by Baekgyu
# Last modified: 200519 by Soomin

SAMPLES = ['HCT116_sDro1', 'HCT116_FL1', 'HCT116_sDro2', 'HCT116_FL2']

MINQUAL = 20
MINQUALPERCENT = 90
ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'
MINLEN = 18
MAXLEN = 26
RANDLEN = 4


rule all:
	input:
		expand('scratch_qc/{sample}.qstats.txt', sample=SAMPLES),
		expand('scratch_qc/{sample}.png', sample=SAMPLES),
		expand('scratch_qc/{sample}.{extension}', sample=SAMPLES, extension=['qfilt.txt', 'fasta.gz']),
		expand('scratch_trim/{sample}.rm3ad.{extension}', sample=SAMPLES, extension=['txt', 'fasta.gz']),
		expand('scratch_trim/{sample}.trimmed.{extension}', sample=SAMPLES, extension=['txt', 'fasta.gz']),
		expand('scratch_spikeins/{sample}.sam', sample=SAMPLES),
		expand('scratch_spikeins/{sample}_SpikeIn_uniq.txt', sample=SAMPLES),
		expand('scratch_map/{sample}.nonspikelist', sample=SAMPLES),
		expand('scratch_map/{sample}_spikefiltered.fasta', sample=SAMPLES),
		expand('scratch_map/{sample}.bam', sample=SAMPLES),
		expand('alignments_hg38/{sample}_sorted.bam', sample=SAMPLES),
		expand('alignments_hg38/{sample}_sorted.bam.bai', sample=SAMPLES),
		expand('annotations/{sample}.intersect', sample=SAMPLES),
		expand('scratch_count/{sample}.readcount', sample=SAMPLES),
		expand('scratch_map/{sample}.mapped.txt', sample=SAMPLES)



rule qc_stats:
	input:
		'sequences/{sample}.fastq.gz'
	output:
		'scratch_qc/{sample}.qstats.txt'
	shell:
		'zcat {input} | fastx_quality_stats -Q 33 -> {output}'

rule qc_figures:
	input:
		'scratch_qc/{sample}.qstats.txt'
	output:
		'scratch_qc/{sample}.png'
	params:
		outputprefix='{sample}'
	shell: 
		'fastq_quality_boxplot_graph.sh -i {input} -o {output} -t {params.outputprefix}'

rule quality_filter:
	input:
		'sequences/{sample}.fastq.gz'
	output: log_qfiltered='scratch_qc/{sample}.qfilt.txt',
			file='scratch_qc/{sample}.fasta.gz'
	shell:
		'zcat {input} |'
		'fastq_quality_filter -Q 33 -q {MINQUAL} -p {MINQUALPERCENT} -v 2>> {output.log_qfiltered} | '
		'fastq_to_fasta -Q 33 -n -r | gzip -c - > {output.file}'

rule remove_3ADseq:
	input: 
		'scratch_qc/{sample}.fasta.gz'
	output: 
		log_3adremoved='scratch_trim/{sample}.rm3ad.txt',
		file='scratch_trim/{sample}.rm3ad.fasta.gz'
	shell: 
		'zcat {input} |'
		'fastx_clipper -Q 33 -a {ADAPTER} -n -c -v 2>> {output.log_3adremoved} |'
		'gzip -c - > {output.file}'

rule trim_random4mer:
	input:
		'scratch_trim/{sample}.rm3ad.fasta.gz'
	output:
		log_trimmed='scratch_trim/{sample}.trimmed.txt',
		file='scratch_trim/{sample}.trimmed.fasta.gz'
	shell:
		'zcat {input} | '
		'cutadapt -u {RANDLEN} -u -{RANDLEN} -m {MINLEN} -M {MAXLEN} - 2>> {output.log_trimmed} | '
		'gzip -c - > {output.file}'


rule align_to_SpikeIns:
	input:
		'scratch_trim/{sample}.trimmed.fasta.gz'
	output:
		'scratch_spikeins/{sample}.sam'
	threads: 16
	shell:
		'STAR '
		'--genomeDir ~/genome/AQseq_spkn/genome.star/ '
		'--readFilesIn {input} --runThreadN {threads} '
		'--readFilesCommand zcat --outSAMunmapped Within '
		'--outFilterMultimapNmax 30 '
		'--outFilterMismatchNmax 1 '
		'--outFileNamePrefix scratch_spikeins/ '
		'--outStd SAM > {output}'

	# For mapping to SpikeIns,
	# '--outSAMunmapped Within: output unmapped reads within the main SAM file'
	# '--outFilterMultimapNmax 30: Allow multiple mapping to 30 SpikeIns'
	# '--outFilterMismatchNmax 1: Allow one mismatch'

rule filter_SpikeIn_uniq_reads:
	input:
		'scratch_spikeins/{sample}.sam'
	output:
		'scratch_spikeins/{sample}_SpikeIn_uniq.txt'
	shell:
		"samtools view -@ 4 -F 4 {input} | grep 'NH:i:1' |"
		"cut -f 3 | sort | uniq -c | tr -d ' ' | tr 'spkn' '\tspkn'> {output}"

	# -F 4: get reads mapped to SpikeIns

rule generate_nonSpikeIn_read_list:
	input:
		'scratch_spikeins/{sample}.sam'
	output:
		'scratch_map/{sample}.nonspikelist'
	threads: 2
	shell:
		'samtools view {input} | cut -f 1 | uniq > {output}'

	# nonspikelist actually contains all reads including spike-ins, but it doesn't matter because spike-ins are designed not to be mapped the human miRNAs.

rule generate_filtered_fastq:
	input:
		nonspikelist='scratch_map/{sample}.nonspikelist',
		original='scratch_trim/{sample}.trimmed.fasta.gz'
	output:
		'scratch_map/{sample}_spikefiltered.fasta'
	shell:
		'faSomeRecords <(zcat {input.original}) {input.nonspikelist} {output}'

rule align_to_human_genome:
	input:
		'scratch_map/{sample}_spikefiltered.fasta'
	output:
		'scratch_map/{sample}.bam'
	threads: 32
	shell: 
		"STAR --genomeDir ~/genome/hg38/gencode.genome.star/ "
		"--readFilesIn {input} --runThreadN {threads} "
		"--readFilesCommand cat --outFilterMultimapNmax 20 "
		"--outFilterMismatchNmax 0"
		"--alignIntronMax 1 "
		"--outFileNamePrefix scratch_map/ "
		"--outStd SAM | samtools view -@ 4 -bS -o {output} -" 

	# For mapping to human genome
	# --outFilterNultimapNmax 20: allow multiple mapping upto 20 loci
	# --outFilterMismatchNmax 0: no mismatch is allowed
	# --alignIntronMax 1: not consider spliced reads 

	# -@ samtools view: multi-threaded # of cores 

rule sort_bam:
	input: 
		'scratch_map/{sample}.bam'
	output:
		'alignments_hg38/{sample}_sorted.bam'
	threads: 32
	shell:
		'samtools sort -@ {threads} {input} -o {output}'

rule index_bam:
	input:
		'alignments_hg38/{sample}_sorted.bam'
	output:
		'alignments_hg38/{sample}_sorted.bam.bai'
	shell:
		'samtools index {input}'

rule intersect_anno:
	input:
		'alignments_hg38/{sample}_sorted.bam'
	output:
		'annotations/{sample}.intersect'
	shell:
		'bedtools intersect -bed -abam {input} '
		'-b ~/genome/hg38/trampoline/annotations.bed.gz '
		'-wa -wb -s -split > {output}'

rule conut_miRNA_reads:
	input:
		'alignments_hg38/{sample}_sorted.bam'
	output:
		'scratch_count/{sample}.readcount'
	shell:
		'bedtools intersect -a ~/genome/hg38/mirbase/mature_pri_BK.bed '
		'-b {input} -wa -s -c > {output}'


rule count_genome_mapped:
	input: 
		'scratch_map/{sample}.bam'
	output:
		'scratch_map/{sample}.mapped.txt'
	shell:
		'samtools view -F 4 {input} | wc -l > {output}'