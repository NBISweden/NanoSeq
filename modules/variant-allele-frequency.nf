process VARIANT_ALLELE_FREQUENCY {

	// Directives

	debug false
	tag "${meta.id}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(mut), path(indel), path(cov), path(duplex_cram)

	output:
	tuple val(meta), path("${meta.id}.vcf.gz*"), emit: ch_vaf_out
	tuple val(task.process), val('bcftools'), eval('bcftools --version | head -n 1 | sed "s/.* //"'), topic: versions

	script:
	"""

		snv_merge_and_vaf_calc.R ${mut[0]} ${indel[0]} ${duplex_cram[0]} ${cov[0]} ${meta.id}.vcf
		bcftools sort -Oz ${meta.id}.vcf -o ${meta.id}.vcf.gz
		bcftools index -t ${meta.id}.vcf.gz
		rm ${meta.id}.vcf

	"""

}
