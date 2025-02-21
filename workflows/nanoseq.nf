/* 
----------------------------------------------------------------------------------------
Main workflow
----------------------------------------------------------------------------------------
*/

// Feature flags

	nextflow.preview.topic = true
	nextflow.preview.output = true

// Import modules

	include { INDEX_REFERENCE } from '../modules/index.nf'
	include { ADD_NANOSEQ_FASTQ_TAGS } from '../modules/add-nanoseq-fastq-tags.nf'
	include { BWA_MEM2_MAP } from '../modules/bwa-map.nf'
	include { REMAP_SPLIT } from '../modules/remap-split.nf'
	include { BWA_MEM2_REMAP } from '../modules/bwa-remap.nf'
	include { MARK_DUPLICATES } from '../modules/mark-duplicates.nf'

// Main workflow

	workflow NANOSEQ {

		take:

			ch_samplesheet
			ch_reference

		main:

			// Branch samplesheet channel according to analysis flow

				ch_samplesheet
					.map { meta, files ->
						return [meta, files]
					}
					.branch { meta, files ->
						fastq: meta.format == 'fastq'
						bam_map: meta.format == 'bam' && meta.mapping // TODO: revisit if channel is needed, or can be combined with cram + renamed with both
						bam_nomap: meta.format == 'bam' && !meta.mapping // TODO: revisit if channel is needed, or can be combined with cram + renamed with both
						cram_map: meta.format == 'cram' && meta.mapping
						cram_nomap: meta.format == 'cram' && !meta.mapping
					}
					.set { ch_samplesheet }

			// Create BWA and samtools index files if they don't exist

				INDEX_REFERENCE (ch_reference)

			// Conditional processing of inputs

				// FASTQ steps

					// Create separate FASTQ channels of duplex and normal reads for parallel processing

						ch_samplesheet.fastq
							.multiMap { meta, files ->
								meta_duplex = meta + [type: "duplex"]
								meta_normal = meta + [type: "normal"]
								duplex: [meta_duplex, [files[0], files[1]]]
								normal: [meta_normal, [files[2], files[3]]]
							}
							.set { ch_fastqs }

					// Add NanoSeq tags to FASTQ, take both duplex and normal read channels

						ADD_NANOSEQ_FASTQ_TAGS (ch_fastqs.duplex.mix(ch_fastqs.normal))

					// Map FASTQ samples

						BWA_MEM2_MAP (ADD_NANOSEQ_FASTQ_TAGS.out.ch_tagged_fastqs, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

				// CRAM steps (NOTE: this includes CRAM output from FASTQ mapping, or user input with/without indexes)
					//FIXME: possibly CRAM and BAM are actually treated the same in the original, but not explicitly stated - revisit and state if so

					// User CRAM input requiring remapping (no indexes) //FIXME: currently not supported, there are some CPAN libraries missing and we either need to refactor process behaviour with Sanger input, or ask them to update the conda libraries

						ch_samplesheet.cram_map
							.multiMap { meta, files ->
								meta_duplex = meta + [type: "duplex"]
								meta_normal = meta + [type: "normal"]
								duplex: [meta_duplex, [files[0]]]
								normal: [meta_normal, [files[1]]]
							}
							.set { ch_cram_remap }

					// CRAM remapping steps

						//REMAP_SPLIT (ch_cram_remap.duplex.mix(ch_cram_remap.normal), ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect()) //FIXME: see above

						//BWA_MEM2_REMAP (REMAP_SPLIT.out.TODO:, ch_reference.collect()) //FIXME: see above

					// User CRAM input, no remapping (with indexes) TODO:

						//ch_samplesheet.cram_nomap
						//	.multiMap { meta, files ->

					// User FASTQ input to CRAM intermediate - mark duplicates //FIXME: later this probably needs to be mixed with user CRAM no index input (.mix(BWA_MEM2_REMAP.out.ch_cram))

						MARK_DUPLICATES (BWA_MEM2_MAP.out.ch_cram, ch_reference)









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
