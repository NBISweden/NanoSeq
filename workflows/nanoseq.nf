/* 
----------------------------------------------------------------------------------------
Main workflow
----------------------------------------------------------------------------------------
*/

// Feature flags

	nextflow.preview.topic = true
	nextflow.preview.output = true

// Import modules

	include { BWA_MEM2_INDEX } from '../modules/bwa-index.nf'
	include { SAMTOOLS_INDEX } from '../modules/samtools-index.nf'
	include { ADD_NANOSEQ_FASTQ_TAGS } from '../modules/add-nanoseq-fastq-tags.nf'
	include { BWA_MEM2_MAP } from '../modules/bwa-map.nf'

// Main workflow

	workflow NANOSEQ {

		take:

			ch_samplesheet
			ch_reference

		main:

			// Create BWA index files if they don't exist

				BWA_MEM2_INDEX(ch_reference)

			// Create samtools index if it doesn't exist

				SAMTOOLS_INDEX(ch_reference)

			// Add NanoSeq tags to FASTQ

				ADD_NANOSEQ_FASTQ_TAGS(ch_samplesheet)

			// Map samples

				BWA_MEM2_MAP(ADD_NANOSEQ_FASTQ_TAGS.out.ch_tagged_fastqs, ch_reference, BWA_MEM2_INDEX.out.ch_bwa_indexes)

			//




			// Report package versions

				Channel.topic('versions')
					.map { process, tool, version ->
						return [process: process, tool: tool, version: version]
					}
					.unique()
					.collect()
					.map { it.join('\n') }
					.collectFile(name: 'package_versions.txt', newLine: true)
					.set { ch_versions }

		emit:

			// Emit channels for publication

				ch_versions = ch_versions

	}
