# Changelog

## [Unreleased]

### Planned

- add support for unindexed bam/cram files

___

<!-- Release history -->


## [1.0.0-NBIS] - TODO: DATE

### Dependencies

- pixi

### Changed or removed

#### Installation & dependencies

- Original repository provided two main installation options: either a large docker image containing all dependencies, or a manual installation option, both have been removed
- Base dependencies are now handled by `pixi` (see `pixi.toml`). Most importantly, `pixi` specifies the Nextflow and Apptainer versions, and fixes their dependencies via the provided `pixi.lock` file
- After cloning the repo, installation is done by:

```
curl -fsSL https://pixi.sh/install.sh | bash

pixi install
```

- Further dependencies within processes are now handled by container images, pulled on the fly by Nextflow & Apptainer
- Where possible, the Seqera containers service was used to build these (i.e., when all dependencies were available via conda channels)
- A gitpod.yml is now provided for small scale testing in a cloud development environment

#### Src

- Removed `src`
- The C/C++ source code is now precompiled in a publically hosted docker image: [GitHub](https://github.com/NBISweden/NanoSeq-src-Docker), [DockerHub](https://hub.docker.com/repository/docker/cormackinsella/nanoseq-src/general), and pulled on the fly

#### Tool versions

- Tool versions have been updated in many cases, including in the `src` Docker image

#### Samplesheet input

- The original workflow asked users to edit the column headers of their samplesheet depending on whether `FASTQ`, `CRAM`, or `BAM` was provided as input. This prevented mixing of input formats, and added potential for user error
- Refactored to a fixed samplesheet layout, and added two columns (`input_format`, `normal_method`), in which users state the file type provided per sample, and the method used for sequencing the matched normal sample (duplex or standard)
- The initialisation subworkflow now loads the samplesheet and validates essential information (unique sample ids, valid input formats, and matched normal sequencing method). If passing checks, file paths are added according to the input type. This is emitted as a channel to the main workflow, where it is branched into individual channels for each input type. The existing logic that allowed parallel mapping for normal and duplex read-pairs has been maintained via a equivalent multimap operation
- Users can now mix input types in a single run

#### Support for unindexed bam/cram files

- Removed support for unindexed bam/cram file inputs as these processes relied on internal Sanger perl libraries and scripts

#### Reference genome & indexes

- Added support for additional FASTA extensions
- Removed hard-coded genome file name, accepts any correctly formatted FASTA file with any prefix, with a rational default
- The original workflow required indexing of the reference genome prior to the workflow
- This has been replaced with a process that executes only if it does not find the indexes, generating them. The process stores the indexes alongside the input FASTA for subsequent runs

#### Output publication

- Outputs are now published via the newer Nextflow output feature

#### Package version reporting

- Package version reporting within process script blocks was replaced with a topic channel, published in the same way as other outputs

#### Overall layout

- Repository layout has been simplified and modularised
- Entry workflow now defined in `main.nf`
- Main workflow now defined in `workflows/nanoseq.nf`
- Initialisation subworkflow defined in `subworkflows/initialise.nf`
- Processes are now one per file and found in `modules`
- Other named workflows have been removed and reimplemented as modules/processes
- Assorted scripts are now found in `bin`
- Basic nextflow config now found in `nextflow.config`
- Parameters are now exclusively defined in `conf/parameters.config`, these have been reorganised and annotated
- Additional profiles have been defined, and all are now found in `conf/profiles.config`
- any `ext.arg` variables are now found in `conf/modules.config`
- Removed redundant configs/deprecated DSL statements

#### Hard coding/system dependency

- Where found, hard coding has been removed and replaced with Nextflow variables or parameters
- Where system paths were used, these have been replaced to be system agnostic





<!-- Pre-fork release history -->

# Releases below are from [fa8sanger/NanoSeq](https://github.com/fa8sanger/NanoSeq), ending at commit 31e34bf

## 3.5.7
* Update to discarded variant script to allow sorting of vcf discarded variant file

## 3.5.6
* Adding parameters to "bcftools mpileup" call to resolve issue for samples using multiple wells
* Correcting code version

## 3.5.5
* Fix bugs in nanoseq result plotter
* Fix bugs in indel caller

## 3.5.4
* Fix bugs in nanoseq result plotter

## 3.5.3
* Fix bugs in nanoseq result plotter

## 3.5.2
* Remove non-deterministic behaviour from indel caller
* Missing semicolon in snv_merge_and_vaf_calc.R

## 3.5.1
* Dockerfile updated to avoid snv_merge_and_vaf_calc.R bug and to improve reproducibility
* Small change to hardcoded overlapping_mask value in indel caller

## 3.5.0
* Minor bug in the indel pipeline fixed (was using bitwise OR instead of logical OR)
* New quality metrics added to the indel calls such that SNVs and indels come out with the same metrics

## 3.4.0

* Fixed faulty PART step that was dropping last genomic intervals for targeted experiments
* Refactored how indels get merged to avoid errors when scaling up calculations
* Updated htslibs, samtools, bcftools and libdeflate

## 3.3.0

* Added new INFO fileds to the final vcf
* Changes to indel calling to reduce false calls
* Added burdens output to nanoseq plotter
* Fixed bugs in the VAF script

## 3.2.1

* Fixed bug in part with empty intervals at end of genome

## 3.2.0

* Modified BAM processing as to keep unmapped reads
* Fixed VCF filter bug in VAF script 

## 3.1.0

* Fixed type error of VAF script

## 3.0.0

* Made changes to accomodate targeted NanoSeq

## 2.3.3

* Added quotes add various places to support ant chromosomes
* Fixed error with 1 bp intervals
* Fixed efficiency calculation problem with coordinates outside of chr ends

## 2.3.2

* Reverted to old efficiency calculation

## 2.3.1

* Added safe creation of tmp dir
* Re-read input file to reduce mem usage in efficiency calculation

## 2.3.0

* Cleaned up/ simplified variantcaller.R
* Added extra check for turncated dsa output

## 2.2.0

* Tested with targeted data
* Fixed issue with coverage done files
* Fixed bamcov so that it handles chromosomes without coverage
* Fixed required/optional arguments from docustring
* Fixed incorrect default argument (-d)
* Added argument for file name of results
* Removed plugin option from htslib that was causing issues with CRAM
* Set seqs_per_slice to 1000 in htslib for improved CRAM streaming 
* Fixed off by one errors in part
* Fixed interval test in part
* Added consistency checks for contigs in BAM files 

## 2.1.0

* Added CRAM support
* Fixed issue with sort in indel pipeline

## 2.0.0

* Fixed 1-off error in dsa code.
* Incorporated Indel scripts.
* Incorporated new code for quick coverage estimation.
* Rewrote Python wrapper.
* Updated licenses.

## 1.5.3

* Fixed bug with last interval for outlier case

## 1.5.1

* Added more error captures for htslib calls

## 1.5.0

* Updated R library installation script
* Added option to skip post processing step
* Simplified excluded chromosome option
* Reduced default coverage window

## 1.4.0

* Initial public release
* Updated to htslib 1.11
* Added licencing information
* Added error checks when reading
