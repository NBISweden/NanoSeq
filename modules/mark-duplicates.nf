process MARK_DUPLICATES {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_low'
	container 'oras://community.wave.seqera.io/library/biobambam_htslib_samtools:3279af2af9e82030'

	// I/O & script

	input:
	tuple val(meta), path(cram)
	path reference_fasta

	output:
	tuple val(meta), path("${meta.id}_${meta.type}.mark.cram*"), emit: ch_cram
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions
	tuple val(task.process), val('biobambam2'), eval('bamsormadup --version 2>&1 | head -n 1 | sed "s/.*version //; s/.$//"'), topic: versions

	script:
	"""

	# Mark duplicates

		samtools view --threads ${task.cpus} --with-header --no-PG ${cram} | \
		bamsormadup inputformat=sam level=0 blocksortverbose=0 rcsupport=1 threads=${task.cpus} fragmergepar=${task.cpus} optminpixeldif=10 | \
		bammarkduplicatesopt verbose=0 level=0 index=0 optminpixeldif=2500 | \
		samtools view --threads ${task.cpus} --no-PG --output ${meta.id}_${meta.type}.mark.cram

	# Index the CRAM

		samtools index --threads ${task.cpus} ${meta.id}_${meta.type}.mark.cram

	"""

}
