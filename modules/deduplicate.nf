process DEDUPLICATE {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_single'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(cram)
	path reference_fasta

	output:
	tuple val(meta), path("${meta.id}_${meta.type}.neat.cram*"), emit: ch_cram
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

	NLINES1=\$(samtools view -H ${cram[0]} | grep ^@PG | grep ID:bamaddreadbundles | wc -l) || true
	NLINES2=\$(samtools view -H ${cram[0]} | grep ^@PG | grep ID:randomreadinbundle | wc -l) || true

	if [ \$NLINES1 -gt 0 ] && [ \$NLINES2 -eq 0 ]

		then
			randomreadinbundle -I ${cram[0]} -O ${meta.id}_${meta.type}.neat.cram -m ${params.minimumReadsInBundle}
			samtools index ${meta.id}_${meta.type}.neat.cram

		else

			cp ${cram[0]} ${meta.id}_${meta.type}.neat.cram
			cp ${cram[1]} ${meta.id}_${meta.type}.neat.cram.crai

	fi

	"""

}
