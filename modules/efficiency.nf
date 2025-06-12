process EFFICIENCY {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_low'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	input:
	tuple val(meta), path(files)
	path reference_fasta
	path indexes

	output:
	path "${meta.id}_${meta.type}.efficiency*tsv", emit: ch_efficiency_tsv
	path "${meta.id}_${meta.type}.efficiency.RBs.pdf", emit: ch_efficiency_pdf

	script:
	"""

		efficiency_nanoseq.pl -t ${task.cpus} -d ${files[2]} -x ${files[0]} -o ${meta.id}_${meta.type}.efficiency -r ${reference_fasta}

	"""

}
