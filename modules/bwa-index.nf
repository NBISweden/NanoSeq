process BWA_MEM2_INDEX {

	// Directives

	debug false
	tag "${reference_fasta}"
	label 'process_low'
	container 'oras://community.wave.seqera.io/library/bwa-mem2:2.2.1--e269358148d98816'
	publishDir path: params.referencePath, mode: 'copy', overwrite: true

	// I/O & script

	input:
	path (reference_fasta)

	output:
	path("${reference_fasta}.*"), emit: ch_bwa_indexes
	tuple val(task.process), val('bwa-mem2'), eval('bwa-mem2 version 2>/dev/null'), topic: versions

	script:
	"""

	bwa-mem2 index ${reference_fasta}

	"""

}
