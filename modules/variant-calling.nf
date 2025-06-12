process VARIANT_CALLING {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}_${meta.jobindex}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams), path(dsa_bed)
	path reference_fasta
	path indexes

	output:
	// tuple  val(metaOut), path(duplex), path(index_duplex), path(normal), path( index_normal), emit : done
	// path "var/${ii}.done"
	// path "var/${ii}.var"
	// path "var/${ii}.cov.bed.gz"
	// path "var/${ii}.discarded_var"
	// path 'var/nfiles' optional true
	// path 'var/args.json' optional true

	script :
	"""

		# Rather than require the dsa nfiles (no longer created), get the value directly from the parameter

			echo ${params.jobs} > nfiles

		# Run variant calling script

			nanoseq.py --ref ${reference_fasta} --duplex ${crams[0]} --normal ${crams[2]} --index ${meta.jobindex} --max_index ${params.jobs} var \
				-a ${params.min_as_xs} \
				-b ${params.min_matched_normal_reads} \
				-c ${params.clips_fraction} \
				-d ${params.min_duplex_depth} \
				-f ${params.consensus_fraction} \
				-i ${params.indel_fraction} \
				-m ${params.min_cycle_number} \
				-n ${params.max_mismatches} \
				-p ${params.proper_pairs_fraction} \
				-q ${params.min_consensus_base_quality} \
				-r ${params.post_trim_read_length} \
				-v ${params.max_normal_vaf} \
				-x ${params.max_cycle_number} \
				-z ${params.min_normal_coverage}

	"""

}
