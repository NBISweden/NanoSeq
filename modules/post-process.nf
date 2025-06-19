process POST_PROCESS {

	// Directives

	debug false
	tag "${meta.id}"
	label 'process_low'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams), path(var1), path(var2), path(var3), path(indel1), path(indel2)
	path reference_fasta
	path indexes

	output:
	tuple val(meta), path("${meta.id}.muts.vcf.gz*"), path("${meta.id}.indel.vcf.gz*"), path("${meta.id}.cov.bed.gz*"), emit: ch_post_process
	path("*.csv"), emit: ch_csv
	path("*.tsv"), emit: ch_tsv, optional: true
	path("*.pdf"), emit: ch_pdf, optional: true

	script:
	def args = task.ext.args ?: ''

	"""

	# Perform logic check that all files have been correctly received

		echo ${params.jobs} > nfiles

		FILES_EXPECTED=\$(cat nfiles || true)
		NUMBER_COV_BED=\$(ls *.cov.bed.gz 2>/dev/null | wc -l || true)
		NUMBER_DISCARDED_VAR=\$(ls *.discarded_var 2>/dev/null | wc -l || true)
		NUMBER_INDEL_FILTERED=\$(ls *.indel.filtered.vcf.gz 2>/dev/null | wc -l || true)
		NUMBER_INDEL_FILTERED_TBI=\$(ls *.indel.filtered.vcf.gz.tbi 2>/dev/null | wc -l || true)
		NUMBER_VAR=\$(ls *.var 2>/dev/null | wc -l || true)

		if [[ "\$FILES_EXPECTED" -ne "\$NUMBER_COV_BED" || "\$FILES_EXPECTED" -ne "\$NUMBER_DISCARDED_VAR" || "\$FILES_EXPECTED" -ne "\$NUMBER_INDEL_FILTERED" || "\$FILES_EXPECTED" -ne "\$NUMBER_INDEL_FILTERED_TBI" || "\$FILES_EXPECTED" -ne "\$NUMBER_VAR" ]]; then
			echo "ERROR: Expected \$FILES_EXPECTED of each file type. Found \$NUMBER_COV_BED cov.bed.gz files, \$NUMBER_DISCARDED_VAR discarded_var files, \$NUMBER_INDEL_FILTERED indel.filtered.vcf.gz files, \$NUMBER_INDEL_FILTERED_TBI indel.filtered.vcf.gz.tbi files, and \$NUMBER_VAR var files."
			exit 1
		fi

	# Run post-processing script

		nanoseq.py --ref ${reference_fasta} --duplex ${crams[0]} --normal ${crams[2]} --threads ${task.cpus} post --name ${meta.id} ${args}

	"""

}
