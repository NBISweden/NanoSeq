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
					.branch {
						fastq: { meta, files -> meta.format == 'fastq' }
						bam_map: { meta, files -> meta.format == 'bam' && meta.mapping }
						bam_nomap: { meta, files -> meta.format == 'bam' && !meta.mapping }
						cram_map: { meta, files -> meta.format == 'cram' && meta.mapping }
						cram_nomap: { meta, files -> meta.format == 'cram' && !meta.mapping }
					}
					.set { ch_samplesheet }

			// Create BWA and samtools index files if they don't exist

				INDEX_REFERENCE(ch_reference)

			// Conditional processing of inputs

				// FASTQ

					// Multimap the FASTQ input channel to emit duplex and normal reads separately

						ch_samplesheet.fastq
							.multiMap { meta, files ->
								meta_duplex = meta.clone() // Clone metadata without editing original
								meta_normal = meta.clone()
								meta_duplex["type"] = "duplex" // Add type to metadata
								meta_normal["type"] = "normal"
								duplex: [meta_duplex, [files[0], files[1]]] // Emit read types as separate channels for parallel processing
								normal: [meta_normal, [files[2], files[3]]]
							}
							.set { ch_fastqs }

					// Add NanoSeq tags to FASTQ, both duplex and normal channels

						ADD_NANOSEQ_FASTQ_TAGS(ch_fastqs.duplex.mix(ch_fastqs.normal))

					// Map FASTQ samples

						BWA_MEM2_MAP(ADD_NANOSEQ_FASTQ_TAGS.out.ch_tagged_fastqs, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())








				// // CRAM

						// TODO: cram_map channel, has two files, no indexes, reads 0 and 1 (duplex - normal)
						//REMAP_SPLIT(ch_samplesheet.cram_map, ch_reference)

					// TAKES: cram input channel, "index_dir"
					// DOES: REMAP_SPLIT on the above
					// REMAP on REMAP_split out + index dir
					// EMITS cram from REMAP
					// INDEX dir, is jsut the path to references
					// ch_cram


				// BAM







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
