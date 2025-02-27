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
	include { MARK_DUPLICATES } from '../modules/mark-duplicates.nf'
	include { ADD_READ_BUNDLES } from '../modules/add-read-bundles.nf'
	include { DEDUPLICATE } from '../modules/deduplicate.nf'
	include { EFFICIENCY } from '../modules/efficiency.nf'

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
						cram_bam: (meta.format == 'cram' || meta.format == 'bam')
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

				ch_branches.cram_bam
					.multiMap {meta, files ->
						def meta_duplex = meta + [type: "duplex"]
						def meta_normal = meta + [type: "normal"]
						duplex: [meta_duplex, [files[0], files[1]]]
						normal: [meta_normal, [files[2], files[3]]]
					}
					.set { ch_cram_bam }

			// Create BWA and samtools index files if they don't exist

				INDEX_REFERENCE (ch_reference)

			// Preprocessing

				// FASTQ only steps

					// Add NanoSeq tags to FASTQ, take both duplex and normal read channels

						ADD_NANOSEQ_FASTQ_TAGS (ch_fastqs.duplex.mix(ch_fastqs.normal))

					// Map FASTQ samples

						BWA_MEM2_MAP (ADD_NANOSEQ_FASTQ_TAGS.out.ch_tagged_fastqs, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

					// Mark duplicates

						MARK_DUPLICATES (BWA_MEM2_MAP.out.ch_cram, ch_reference.collect())

				// Shared steps

					// Add read bundles

						ADD_READ_BUNDLES (MARK_DUPLICATES.out.ch_cram.mix(ch_cram_bam.duplex, ch_cram_bam.normal), ch_reference.collect())

					// Deduplicate

						DEDUPLICATE (ADD_READ_BUNDLES.out.ch_cram, ch_reference.collect())


					// TODO: Verify BAM ID

			// Efficiency calculations: note at this stage channels are joined on experiment and type (i.e. the metadata field)

				// Join, return channel with meta + files in one list

					ADD_READ_BUNDLES.out.ch_cram.join(DEDUPLICATE.out.ch_cram).map { meta, read_bundles, deduplicate -> [ meta, read_bundles + deduplicate ] }.set { ch_normal_efficiency }

				// Run calculations

					EFFICIENCY (ch_normal_efficiency, ch_reference.collect())





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
