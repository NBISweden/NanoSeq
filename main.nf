#!/usr/bin/env nextflow

/* 
----------------------------------------------------------------------------------------

NBISweden/NanoSeq

GitHub: https://github.com/NBISweden/NanoSeq

This workflow is a modified version of the original NanoSeq workflow (https://github.com/cancerit/NanoSeq), and was started from a fork off the original project (https://github.com/fa8sanger/NanoSeq, commit 31e34bf).

Contributors to this version:

- Cormac Kinsella (cormac.kinsella@nbis.se)

The original software was licensed under the GNU Affero General Public License v3.0, as is this modified version.

Changes made to the original codebase are summarised in the CHANGES.md file, and detailed in the repo commit history.

----------------------------------------------------------------------------------------
*/

// Imports

	include { INITIALISE } from './subworkflows/initialise.nf'
	include { NANOSEQ } from './workflows/nanoseq.nf'

// Entry workflow

	workflow {

		main:

			INITIALISE ()
			NANOSEQ (
				INITIALISE.out.ch_samplesheet,
				INITIALISE.out.ch_reference,
			)

		publish:

			efficiency = NANOSEQ.out.ch_efficiency_out
			results_one = NANOSEQ.out.results_one
			results_two = NANOSEQ.out.results_two
			results_three = NANOSEQ.out.results_three
			package_versions = NANOSEQ.out.ch_versions

	}

// Publish outputs

	output {

		efficiency {
			path { meta, one, two, three -> "${meta.id}/efficiency_reports" }
			mode 'copy'
			overwrite false
		}

		results_one {
			path { meta, one, two, three -> "${meta.id}/variant_calling" }
			mode 'copy'
			overwrite false
		}

		results_two {
			path { meta, one, two, three, four, five -> "${meta.id}/variant_calling_summary_files" }
			mode 'copy'
			overwrite false
		}

		results_three {
			path { meta, one -> "${meta.id}/variant_allele_frequency" }
			mode 'copy'
			overwrite false
		}

		package_versions {
			path 'package_versions'
			mode 'copy'
			overwrite false
		}

	}

// Email report

	if (params.email_report) {

		workflow.onComplete {

			// Prepare email content
			def workflow_status = workflow.success ? 'COMPLETED' : 'FAILED'
			def email_address = params.email
			def subject = "NanoSeq workflow run ${workflow.runName}: ${workflow_status}"
			def msg = """
			Pipeline execution summary
			---------------------------
			Run Name     : ${workflow.runName}
			Completed at : ${workflow.complete}
			Duration     : ${workflow.duration}
			Success      : ${workflow.success}
			Exit status  : ${workflow.exitStatus}
			Error report : ${workflow.errorReport ?: 'No errors'}
			"""
			.stripIndent()

			// Send the email
			sendMail(
				to: email_address,
				subject: subject,
				body: msg
			)

		}
	}
