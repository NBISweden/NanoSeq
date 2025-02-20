process INDEX_REFERENCE {

	// Directives

	debug false
	tag "${reference_fasta}"
	label 'process_low'
	container 'oras://community.wave.seqera.io/library/bwa-mem2_htslib_samtools:ada17f4d0d757cc3'
	storeDir { reference_fasta.toRealPath().parent }

	// I/O & script

	input:
	path reference_fasta

	output:
	path("${reference_fasta}.*"), emit: ch_indexes
	// PLANNED: enable topic channel once Nextflow fixes incompatibility with topic channels
	//tuple val(task.process), val('bwa-mem2'), eval('bwa-mem2 version 2>/dev/null'), topic: versions
	//tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

	bwa-mem2 index ${reference_fasta}
	samtools faidx ${reference_fasta}
	samtools dict ${reference_fasta} > ${reference_fasta}.dict

	"""

}
