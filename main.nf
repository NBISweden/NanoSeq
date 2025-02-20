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

			NANOSEQ.out.ch_versions >> 'package_versions'

	}

// Publish outputs

	output {

		package_versions {
			path 'package_versions'
			mode 'copy'
			overwrite false
		}

	}
