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
	"""

	bwa-mem2 mem -C -t ${task.cpus} ${reference_fasta} ${reads} | \
		samtools sort --threads ${task.cpus} -n -O bam -l 0 -m 2G - | \
		samtools view --threads ${task.cpus} --reference ${reference_fasta} --output ${meta.id}_${meta.type}.cram

	"""

}
