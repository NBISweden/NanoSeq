process COVERAGE {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams)
	path reference_fasta
	path indexes

	output:
	tuple path("args.json"), path("cov.bed.gz"), path("gIntervals.dat"), emit: ch_coverage

	script:
	"""

	# Run NanoSeq coverage script

		coverage.py --ref ${reference_fasta} --normal ${crams[2]} --duplex ${crams[0]} --threads ${task.cpus} cov -Q ${params.minimum_duplex_mapq} --exclude \"${params.contig_exclude}\" --include \"${params.contig_include}\" --larger ${params.contig_larger_than}

	"""

}
