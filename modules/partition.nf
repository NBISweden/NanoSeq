process PARTITION {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams), path(cov_args), path(cov_beds), path(g_intervals), path(nfiles)
	path reference_fasta
	path indexes

	output :
	tuple val(meta), path(crams), path("part_args.json"), path("intervalsPerCPU.dat"), emit: ch_partition
	tuple val(task.process), val('python'), eval('python --version | sed "s/.* //"'), topic: versions
	tuple val(task.process), val('nanoseq.py'), eval('nanoseq.py -v'), topic: versions

	script :
	def args = task.ext.args ?: ''
	def args2 = task.ext.args2 ?: ''

	"""

	# Run NanoSeq partition script

		nanoseq.py  --ref ${reference_fasta} --duplex ${crams[0]} --normal ${crams[2]} --threads ${task.cpus} part --jobs ${params.jobs} ${args} ${args2}

	"""

}
