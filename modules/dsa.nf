process DSA {

	// Directives

	debug true
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
	// 	tuple val(metaOut), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
	// 	path "dsa/${ii}.done"
	// 	path "dsa/${ii}.dsa.bed.gz"
	// 	path "dsa/nfiles" optional true
	// 	path "dsa/args.json" optional true
	// tuple val(task.process), val('python'), eval('python --version | sed "s/.* //"'), topic: versions
	// tuple val(task.process), val('nanoseq.py'), eval('nanoseq.py -v'), topic: versions

	script:
	def args = task.ext.args ?: ''
	def args2 = task.ext.args2 ?: ''
	// Add job index number to output metadata for rejoining
	metaOut = meta.clone()
	metaOut["jobindex"] = jobindex

	"""

		nanoseq.py --ref ${reference_fasta} --normal ${crams[2]} --duplex ${crams[0]} --index ${jobindex} --max_index ${params.jobs} dsa -d ${params.minimum_duplex_depth} -q ${params.minimum_normal_base_quality} ${args} ${args2}

	"""

}
