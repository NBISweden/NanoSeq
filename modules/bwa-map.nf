process BWA_MEM2_MAP {

	// Directives

	debug false
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	container 'oras://community.wave.seqera.io/library/bwa-mem2_htslib_samtools:ada17f4d0d757cc3'

	// I/O & script

	input:
	tuple val(meta), path(reads)
	path reference_fasta
	path indexes

	output:
	tuple val(meta), path("${meta.id}_${meta.type}.cram"), emit: ch_cram
	tuple val(task.process), val('bwa-mem2'), eval('bwa-mem2 version 2>/dev/null'), topic: versions
	tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions
	tuple val(task.process), val('htslib'), eval('bgzip --version | head -n 1 | sed "s/.* //"'), topic: versions

	script:
	def args = task.ext.args ?: '-C'
	def args2 = task.ext.args2
	"""

	bwa-mem2 mem ${args} -t ${task.cpus} ${reference_fasta} ${reads} | \
		samtools sort --threads ${task.cpus} ${args2} - | \
		samtools view --threads ${task.cpus} -T ${reference_fasta} -o ${meta.id}_${meta.type}.cram

	"""

}
