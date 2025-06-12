process DSA {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}_${jobindex}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams), path(part_args), path(intervalsPerCPU)
	path reference_fasta
	path indexes
	each jobindex

	output:
	tuple val(metaOut), path(crams), path("*.dsa.bed.gz"), emit: ch_dsa
	tuple val(task.process), val('python'), eval('python --version | sed "s/.* //"'), topic: versions
	tuple val(task.process), val('nanoseq.py'), eval('nanoseq.py -v'), topic: versions

	script:
	def args = task.ext.args ?: ''
	def args2 = task.ext.args2 ?: ''
	// Add job index number to output metadata for rejoining
	metaOut = meta.clone()
	metaOut["jobindex"] = jobindex

	"""

		# Rather than require the cov_args.json file, get the value directly from the parameter

			echo "{\\"Q\\": ${params.min_duplex_mapq}}" > minimum_duplex_mapq.json

		# Run NanoSeq DSA script

			nanoseq.py --ref ${reference_fasta} --duplex ${crams[0]} --normal ${crams[2]} --index ${jobindex} --max_index ${params.jobs} dsa -d ${params.min_duplex_depth} -q ${params.min_normal_base_quality} ${args} ${args2}

	"""

}
