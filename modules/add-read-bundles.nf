process ADD_READ_BUNDLES {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_single'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(cram)
	path reference_fasta

	output:
	tuple val(meta), path("${meta.id}_${meta.type}.rb.cram*"), emit: ch_cram
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

		bamaddreadbundles -I *.cram -O ${meta.id}_${meta.type}.rb.cram
		samtools index ${meta.id}_${meta.type}.rb.cram

	"""

}
