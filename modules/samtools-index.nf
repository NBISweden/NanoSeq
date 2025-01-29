process SAMTOOLS_INDEX {

	// Directives

	debug false
	tag "${reference_fasta}"
	label 'process_low'
	container 'oras://community.wave.seqera.io/library/htslib_samtools:1.21--7c0846afa354d6db'
	publishDir path: params.referencePath, mode: 'copy', overwrite: true

	// I/O & script

	input:
	path (reference_fasta)

	output:
	path("${reference_fasta}.dict"), emit: ch_samtools_index
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

	samtools dict ${reference_fasta} > ${reference_fasta}.dict

	"""

}
