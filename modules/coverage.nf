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
	tuple path("cov_args.json"), path("*cov.bed.gz"), path("gIntervals.dat"), path("nfiles"), emit: ch_coverage

	script:
	"""

	# Run NanoSeq coverage script

		nanoseq.py --ref ${reference_fasta} --normal ${crams[2]} --duplex ${crams[0]} --threads ${task.cpus} cov -Q ${params.minimum_duplex_mapq} --exclude \"${params.contig_exclude}\" --include \"${params.contig_include}\" --larger ${params.contig_larger_than}

	# Ensure all expected outputs were created

	BEDS_EXPECTED=\$(cat nfiles || true)
	NUMBER_BED_FILES=\$(ls *.cov.bed.gz 2>/dev/null | wc -l || true)
	NUMBER_PROCESSED=\$(ls *.coverage_processed 2>/dev/null | wc -l || true)

	if [[ "\$BEDS_EXPECTED" -ne "\$NUMBER_BED_FILES" || "\$BEDS_EXPECTED" -ne "\$NUMBER_PROCESSED" ]]; then
		echo "ERROR: Expected \$BEDS_EXPECTED bed files, found \$NUMBER_BED_FILES, of which \$NUMBER_PROCESSED were successfully processed"
		exit 1
	fi

	"""

}
