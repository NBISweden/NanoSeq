# Dockerfile for Nanoseq custom code

- For maintainability and readability, this Dockerfile has been limited to strict dependencies of NanoSeq custom code & external tools used alongside custom code. Where possible it eliminates the need for external build scripts found in the original project.

- The image is available on the [GitHub Container Registry](https://github.com/NBISweden/NanoSeq/pkgs/container/nanoseq-src), and is pulled automatically by Apptainer during pipeline execution

## Changelog


### 0.93

- Multi stage build to reduce image size and resolve compiler issues

- NOTE: these packages were not yet added

	- bcftools
	- verifybamid

### 0.92

- Include R & dependencies
- Changed gcc version for R package compilation

### 0.91

- Test perl & dependencies

### 0.9

- Overhauled test image, revised up to and including compiled c/cpp code, removing external build scripts

### 0.1

- Stripped down version of the original Dockerfile/Makefile/build scripts, focusing on getting the C/C++ code to run
