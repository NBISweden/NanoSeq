process ADD_READ_BUNDLES {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_single'
	container 'docker://cormackinsella/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(cram)
	path reference_fasta

	output:
	tuple val(meta), path("${meta.id}_${meta.type}.rb.cram*"), emit: ch_cram
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

	NLINES=\$(samtools view ${cram[0]} | head -1 | grep rb: | grep rc: | grep mb: | grep mc: | wc -l) || true

	if [ \$NLINES -ne 0 ]

		then

			bamaddreadbundles -I *.cram -O ${meta.id}_${meta.type}.rb.cram
			samtools index ${meta.id}_${meta.type}.rb.cram

		else

			cp ${cram[0]}  ${meta.id}_${meta.type}.rb.cram
			cp ${cram[1]}  ${meta.id}_${meta.type}.rb.cram.crai

	fi

	"""

}
