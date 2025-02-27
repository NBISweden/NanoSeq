process EFFICIENCY {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_low'
	container 'oras://community.wave.seqera.io/library/perl:5.32.1--9e3c43247be68b3b'

	input:
	tuple val(meta), path(files)
	path reference_fasta

	output:
	path "${meta.id}_${meta.type}.efficiency.tsv", emit: ch_efficiency_tsv

	script:
	"""

	efficiency_nanoseq.pl -t ${task.cpus} -d ${files[2]} -x ${files[0]} -o ${meta.id}_${meta.type}.efficiency -r ${reference_fasta}

	"""

}
