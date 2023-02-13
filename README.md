
# RNA Sequence Kit

This tool provides several functions for analyzing and modifying RNA sequences. The available functions are:

-   `gc_content`: calculates the GC content of an RNA sequence.
-   `codon_usage`: calculates the codon usage of an RNA sequence.
-   `splice_site_prediction`: predicts the splice sites in an RNA sequence.
-   `multiple_sequence_alignment`: performs multiple sequence alignment on a list of RNA sequences using ClustalW, MUSCLE or MAFFT.
-   `remove_introns`: removes introns from an RNA sequence.
-   `add_mirna_binding_sites`: adds specific miRNA binding sites to an RNA sequence.
-   `convert_sequence_format`: converts an RNA sequence from one format to another (either string or FASTA).
-   `export_sequence`: exports an RNA sequence to a file in a specified format.

#### The following functions are analysis functions:

-   `gc_content`
-   `codon_usage`
-   `splice_site_prediction`

####  The following functions are modification functions:

-   `remove_introns`
-   `add_mirna_binding_sites`
-   `convert_sequence_format`

## Requirements

This tool requires the following Python packages:

-   argparse
-   os
-   re
-   Biopython

## Usage

    python rna_tool.py [-h] [-a] [-m] [-e EXPORT] sequence

The `sequence` argument is required and should be the RNA sequence to be processed.

The following optional arguments are available:

-   `-a`, `--analysis`: Perform analysis on the sequence.
-   `-m`, `--modification`: Perform modification on the sequence.
-   `-e EXPORT`, `--export EXPORT`: Export the processed sequence to a file.

## Contribution

If you want to contribute to this tool, you can open an issue or a pull request. All contributions are welcome!

## License

This tool is licensed under the MIT license. 
