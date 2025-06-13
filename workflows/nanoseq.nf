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
	include { COVERAGE } from '../modules/coverage.nf'
	include { PARTITION } from '../modules/partition.nf'
	include { DSA } from '../modules/dsa.nf'
	include { VARIANT_CALLING } from '../modules/variant-calling.nf'
	include { INDEL } from '../modules/indel.nf'
	include { POST_PROCESS } from '../modules/post-process.nf'

// Main workflow

	workflow NANOSEQ {

		take:

			ch_samplesheet
			ch_reference

		main:

		// Channel preparation

			// Create separate channels for each forward/reverse readset (1 duplex pair + 1 matched normal pair), for parallel processing

				ch_samplesheet
					.multiMap { meta, files ->
						def meta_duplex = meta + [type: "duplex"]
						def meta_normal = meta + [type: "normal"]
						duplex: [meta_duplex, [files[0], files[1]]]
						normal: [meta_normal, [files[2], files[3]]]
					}
					.set { ch_fastqs }

		// Preprocessing modules

			// Create BWA and samtools index files if they don't exist, store alongside reference FASTA

				INDEX_REFERENCE (ch_reference)

			// Extracts 3 nucleotide barcodes from R1 and matched R2 reads, tags FASTQ header with rb and mb (read barcode and mate barcode), strips remaining spacer to leave only non-adapter sequence

				ADD_NANOSEQ_FASTQ_TAGS (ch_fastqs.duplex.mix(ch_fastqs.normal))

			// Map reads to reference genome, output CRAM

				BWA_MEM2_MAP (ADD_NANOSEQ_FASTQ_TAGS.out.ch_tagged_fastqs, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

			// Mark duplicate reads, add biobambam rc and mc tags (unclipped coordinate tags)

				MARK_DUPLICATES (BWA_MEM2_MAP.out.ch_cram, ch_reference.collect())

			// Add read bundle tags, remove reads failing QC

				ADD_READ_BUNDLES (MARK_DUPLICATES.out.ch_cram, ch_reference.collect())

			// Deduplicate read bundles, one read per bundle per strand

				DEDUPLICATE (ADD_READ_BUNDLES.out.ch_cram, ch_reference.collect())

			// For efficiency calculation join pre & post deduplication channels on metadata (i.e. sample id & meta_type)

				ADD_READ_BUNDLES.out.ch_cram
					.join(DEDUPLICATE.out.ch_cram)
					.map { meta, read_bundles, deduplicated -> [ meta, read_bundles + deduplicated ] }
					.set { ch_efficiency }

			// Produce efficiency reports

				EFFICIENCY (ch_efficiency, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

		// Main NanoSeq analysis modules

			// For input, join non-deduplicated duplex and deduplicated normal CRAM files per sample

				ADD_READ_BUNDLES.out.ch_cram
					.filter { meta, read_bundles -> meta.type == 'duplex' }
					.map { meta, read_bundles ->
						[ meta.id, meta, read_bundles ]
					}
					.join(
				DEDUPLICATE.out.ch_cram
					.filter { meta, deduplicated -> meta.type == 'normal' }
					.map { meta, deduplicated ->
						[ meta.id, meta, deduplicated ]
					}
					)
					.map { id, duplex_meta, read_bundles, normal_meta, deduplicated ->
						def meta_join = duplex_meta.clone()
						meta_join.type = 'pair'
						[meta_join, read_bundles + deduplicated]
						}
					.set { ch_nanoseq_crams }

			// Coverage

				COVERAGE (ch_nanoseq_crams, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

			// Partition

				PARTITION (COVERAGE.out.ch_coverage, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

			// Create indexes to split work over n jobs

				jobIndexes = Channel.of(1..params.jobs)

			// DSA

				DSA (PARTITION.out.ch_partition, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect(), jobIndexes)

			// Variant calling

				VARIANT_CALLING (DSA.out.ch_dsa, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

			// Indel calling

				INDEL (DSA.out.ch_dsa, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

			// Stage variant and indel calling results together per sample

				VARIANT_CALLING.out.ch_variant_calling
					// Join channels on their metadata (sample, normalMethod & jobindex can vary)
					.join(INDEL.out.ch_indel_calling)
					.map { meta, crams, var1, var2, var3, indel1, indel2 ->
						// Jobindex no longer required as chunks are to be grouped and this will prevent it
						def metaNew = meta.clone()
						metaNew.remove('jobindex')
						[metaNew, crams, var1, var2, var3, indel1, indel2]
					}
					// Group by variable metadata (sample & normalMethod), when each sample finishes all constituent jobs
					.groupTuple(size: params.jobs, sort: 'hash') // N = jobs due to join operation above (i.e. not jobs * 2), hash sort to ensure deterministic order for resume
					.map { meta, cramsList, var1, var2, var3, indel1, indel2 ->
						// Take only the first set of CRAMS from the list per sample (avoid redundant file staging)
						def crams = cramsList[0]
						[meta, crams, var1, var2, var3, indel1, indel2]
					}
					.set { ch_variant_indel }

			// Post-processing

				POST_PROCESS (ch_variant_indel, ch_reference.collect(), INDEX_REFERENCE.out.ch_indexes.collect())

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

				ch_efficiency_tsv = EFFICIENCY.out.ch_efficiency_tsv
				ch_efficiency_pdf = EFFICIENCY.out.ch_efficiency_pdf
				ch_versions = ch_versions

	}
