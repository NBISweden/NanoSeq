process PARTITION {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams)
	path reference_fasta
	path indexes
	tuple path(cov_args), path(cov_bed), path(g_intervals)

	output :
		// tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
		// path 'part/1.done'
		// path 'part/intervalsPerCPU.dat'
		// path 'part/args.json'


	script :
	def args = task.ext.args ?: ''
	def args2 = task.ext.args2 ?: ''

	"""

	# Run NanoSeq partition script

		# partition.py --threads ${task.cpus} --jobs ${params.jobs} ${args} ${args2}

	"""

}
