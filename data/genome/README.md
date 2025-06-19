# Default location for reference genome FASTA and indexes

## Reference file

- Default location and name is: `${projectDir}/data/genome/genome.fa`
- To override this, use the `--fasta` parameter
- Any FASTA file with the correct format (i.e. '.fa', '.fna', '.fasta') can be supplied

## Index files

- The workflow will look for BWA-MEM2 and samtools dict indexes in the same directory as the reference genome
- There is no user requirement to provide these, as they will be generated on the fly if missing and stored alongisde the FASTA for subsequent runs
