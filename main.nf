#!/usr/bin/env nextflow

/* 
----------------------------------------------------------------------------------------

Main workflow

NanoSeq-Fork

GitHub: https://github.com/NBISweden/NanoSeq-fork

Contributors to the modified version:
- Cormac Kinsella (cormac.kinsella@nbis.se)

This workflow is a modified version of the NanoSeq workflow: https://github.com/cancerit/NanoSeq, and was forked from https://github.com/fa8sanger/NanoSeq, commit 31e34bf.

The original software was licensed under the GNU Affero General Public License v3.0, as is this modified version.

Changes made to the original codebase are summarised in the CHANGES.md file, and detailed in the repo commit history.

----------------------------------------------------------------------------------------
*/

include { VERIFY } from './subworkflows/verify.nf'
include { NANOSEQ_FORK } from './workflows/nanoseq-fork.nf'

workflow {

	VERIFY ()
	NANOSEQ_FORK ()

}
