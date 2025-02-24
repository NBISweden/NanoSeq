process BWA_MEM2_REMAP {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	//FIXME: container ''

	// I/O & script

	input:
	tuple val(meta), path(reads)
	path reference_fasta

	output:
	tuple val(meta), path("sort/${meta.name}.cram"), path("sort/${meta.name}.cram.crai"), emit: cram //TODO: remove this crai
	path  "versions.yml", emit: versions

	script:
	def args = task.ext.args ?: '-n -tags BC,QT,mb,rb -b \'-T 30 -Y\''
	def args2 = task.ext.args ?: '-n'
	"""

	mkdir -p sort

	bwa_mem.pl -p bwamem ${args} -bwamem2 -cram -t ${task.cpus} -mt ${task.cpus} -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} ${reads}
	
	#need the mark process so the final cram file gets placed in the correct location
	bwa_mem.pl -p mark -n ${args} -bwamem2  -cram -t ${task.cpus} -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} ${reads}

	samtools sort -@ ${task.cpus} ${args2} -O cram -m 2G -o ./sort/${meta.name}.cram ./bwamem/${meta.name}.cram
	touch sort/${meta.name}.cram.crai

	"""

}
