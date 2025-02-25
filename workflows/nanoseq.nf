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
	include { ADD_READ_BUNDLES } from '../modules/add-read-bundles.nf'
	include { DEDUPLICATE } from '../modules/deduplicate.nf'

// Main workflow

	workflow NANOSEQ {

		take:

			ch_samplesheet
			ch_reference

		main:

			// Branch samplesheet channel according to analysis flow

				ch_samplesheet
					.branch { meta, files ->
						fastq: meta.format == 'fastq'
						cram_bam_nomap: (meta.format == 'cram' || meta.format == 'bam') && !meta.mapping
						cram_bam_map: (meta.format == 'cram' || meta.format == 'bam') && meta.mapping
					}
					.set { ch_branches }

			// Multimap each branch, create separate channels for duplex and normal reads for parallel processing

				ch_branches.fastq
					.multiMap { meta, files ->
						def meta_duplex = meta + [type: "duplex"]
						def meta_normal = meta + [type: "normal"]
						duplex: [meta_duplex, [files[0], files[1]]]
						normal: [meta_normal, [files[2], files[3]]]
					}
					.set { ch_fastqs }

				ch_branches.cram_bam_nomap
					.multiMap {meta, files ->
						def meta_duplex = meta + [type: "duplex"]
						def meta_normal = meta + [type: "normal"]
						duplex: [meta_duplex, [files[0], files[1]]]
						normal: [meta_normal, [files[2], files[3]]]
					}
					.set { ch_cram_bam_no_map }

				// ch_branches.cram_bam_map
				// 	.multiMap { meta, files ->
				// 		def meta_duplex = meta + [type: "duplex"]
				// 		def meta_normal = meta + [type: "normal"]
				// 		duplex: [meta_duplex, [files[0]]]
				// 		normal: [meta_normal, [files[1]]]
				// 	}
				// 	.set { ch_cram_bam_map }

			// Create BWA and samtools index files if they don't exist

				INDEX_REFERENCE (ch_reference)

			// Conditional processing of inputs

				// FASTQ steps

					// Add NanoSeq tags to FASTQ, take both duplex and normal read channels

						ADD_NANOSEQ_FASTQ_TAGS (ch_fastqs.duplex.mix(ch_fastqs.normal))

					// Map FASTQ samples

						BWA_MEM2_MAP (ADD_NANOSEQ_FASTQ_TAGS.out.ch_tagged_fastqs, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

					// User FASTQ input to CRAM intermediate - mark duplicates

						MARK_DUPLICATES (BWA_MEM2_MAP.out.ch_cram, ch_reference.collect())

				// CRAM and BAM steps (NOTE: this includes CRAM output from FASTQ mapping, and CRAM/BAM user input with or without indexes)

					// User CRAM/BAM input requiring remapping (no indexes) //FIXME: currently not supported, there are some CPAN libraries missing and we either need to refactor process behaviour with Sanger input, or ask them to update the conda libraries
					// CRAM/BAM remapping steps

						//REMAP_SPLIT (ch_cram_bam_map.duplex.mix(ch_cram_bam_map.normal), ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect()) //FIXME: see above

						//BWA_MEM2_REMAP (REMAP_SPLIT.out.TODO:, ch_reference.collect()) //FIXME: see above

					// Add read bundles, CRAM/BAM input with indexes joins here, same as FASTQ input MARK_DUPLICATES.out.ch_cram

						ADD_READ_BUNDLES (MARK_DUPLICATES.out.ch_cram.mix(ch_cram_bam_no_map.duplex, ch_cram_bam_no_map.normal), ch_reference.collect())

					// Deduplicate

						//DEDUPLICATE

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
