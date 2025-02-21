process REMAP_SPLIT {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	//container 'oras://community.wave.seqera.io/library/perl-autodie_perl-const-fast_perl-ipc-system-simple_perl:ab8ae01f6fde6abd' // status with basic perl modules working
	container 'oras://community.wave.seqera.io/library/perl-autodie_perl-capture-tiny_perl-const-fast_perl-ipc-system-simple_pruned:57e190caedddd6b3' // FIXME: now trying to get PCAP libraries
	//container 'oras://community.wave.seqera.io/library/bwa-mem2_perl-autodie_perl-ipc-system-simple_perl:5a1c4e7f74efcaec' FIXME: once basic perl modules working, replace PCAP modules with conda

	// I/O & script

	input:
	tuple val(meta), path(reads)
	path reference_fasta
	path indexes

	output:
	tuple val(meta), path("chunks"), path(cram), emit: ch_cram_remap_chunks // TODO: third element seems superfluous but unsure. Second might be improved too.

	script:
	def args = task.ext.args ?: '-tags BC,QT,mb,rb -b \'-T 30 -Y\''
	"""

	mkdir -p chunks
	bin-bwa_mem.pl -p setup -threads ${task.cpus} ${args} -nomarkdup -bwamem2 -outdir ./chunks -reference ${reference_fasta} -sample ${meta.id}_${meta.type} ${reads}
	bin-bwa_mem.pl -p split -threads ${task.cpus} ${args} -nomarkdup -bwamem2 -cram -outdir ./chunks -reference ${reference_fasta} -sample ${meta.id}_${meta.type} ${reads}

	"""

}
