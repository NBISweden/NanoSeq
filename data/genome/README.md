# Default location for reference genome FASTA and indexes

## Reference file

- Unless specified with `--referencePath`, this location will be searched for the reference assembly file
- Unless specified with `--fasta`, the workflow will expect the file to be named `genome.fa` 
- Any FASTA file with the correct format (i.e. '.fa', '.fna', '.fasta') can be supplied to `--fasta`

## Index files

- The workflow will look for BWA-MEM2 and samtools dict indexes in the same directory as the reference genome, i.e. `--referencePath`
- There is no user requirement to provide these, as they will be generated on the fly if missing, and also copied to `--referencePath` for subsequent runs
