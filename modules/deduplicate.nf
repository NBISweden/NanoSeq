process DEDUPLICATE {

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
	tuple val(meta), path("${meta.id}_${meta.type}.dedup.cram*"), emit: ch_cram
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

		randomreadinbundle -I ${cram[0]} -O ${meta.id}_${meta.type}.dedup.cram -m ${params.minReadsInBundle}
		samtools index ${meta.id}_${meta.type}.dedup.cram

	"""

}
