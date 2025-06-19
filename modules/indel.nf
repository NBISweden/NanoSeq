process INDEL {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}_${meta.jobindex}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams), path(dsa_bed)
	path reference_fasta
	path indexes

	output:
	tuple val(meta), path("${meta.jobindex}.indel.filtered.vcf.gz"), path("${meta.jobindex}.indel.filtered.vcf.gz.tbi"), emit: ch_indel_calling
	tuple val(task.process), val('python'), eval('python --version | sed "s/.* //"'), topic: versions
	tuple val(task.process), val('nanoseq.py'), eval('nanoseq.py -v'), topic: versions

	script:
	"""

	# Rather than require the dsa nfiles (no longer created), get the value from the parameter

		echo ${params.jobs} > nfiles

	# Run indel script

		nanoseq.py --ref ${reference_fasta} --duplex ${crams[0]} --normal ${crams[2]} --index ${meta.jobindex} --max_index ${params.jobs} indel --rb ${params.indel_min_bundle_reads} --t3 ${params.trim_3} --t5 ${params.trim_5} -z ${params.min_normal_coverage} -v ${params.max_normal_vaf} -a ${params.min_as_xs} -c ${params.clips_fraction}

	"""

}
