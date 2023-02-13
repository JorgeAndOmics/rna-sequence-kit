import argparse
import os
import re
from Bio.Align.Applications import ClustalwCommandline

def gc_content(sequence):
    """
    Calculates the GC content of a RNA sequence.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be analyzed.

    Returns
    -------
    float
        The GC content of the sequence.
    """
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = gc_count / len(sequence)
    return gc_content

def codon_usage(sequence):
    """
    Calculates the codon usage of a RNA sequence.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be analyzed.

    Returns
    -------
    dict
        A dictionary where the keys are codons and the values are their corresponding usage count in the sequence.
    """
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    codon_usage = {}
    for codon in codons:
        if codon in codon_usage:
            codon_usage[codon] += 1
        else:
            codon_usage[codon] = 1
    return codon_usage

def splice_site_prediction(sequence):
    """
    Predicts splice sites in a RNA sequence.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be analyzed.

    Returns
    -------
    list
        A list of the predicted splice sites in the sequence.
    """
    splice_site_pattern = re.compile(r'(GT|AG)')
    return [m.start() for m in splice_site_pattern.finditer(sequence)]

def multiple_sequence_alignment(sequences, algorithm='ClustalW'):
    """
    Performs multiple sequence alignment on a list of RNA sequences.

    Parameters
    ----------
    sequences : list
        A list of RNA sequences to be aligned.
    algorithm : str, optional
        The algorithm to be used for the multiple sequence alignment, either 'ClustalW', 'MUSCLE', or 'MAFFT' (default is 'ClustalW').

    Returns
    -------
    list
        A list of the aligned sequences.
    """
    if algorithm == 'ClustalW':
        from Bio.Alphabet import generic_dna
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO

        records = [SeqRecord(seq, id=str(i)) for i, seq in enumerate(sequences)]
        input_file = "input.fa"
        output_file = "output.aln"
        SeqIO.write(records, input_file, "fasta")

        cmd = ClustalwCommandline("clustalw2", infile=input_file, outfile=output_file)
        cmd()

        alignment = list(SeqIO.parse(output_file, "clustal"))
        return alignment

    else:
        raise ValueError("Invalid algorithm. Choose either 'ClustalW', 'MUSCLE', or 'MAFFT'.")

def remove_introns(sequence):
    """
    Removes introns from a RNA sequence.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be modified.

    Returns
    -------
    str
        The modified RNA sequence with introns removed.
    """
    exons = re.findall(r'(AUG[A-Z]{3})*?(?=(AUG))', sequence)
    return "".join(exons)

def add_mirna_binding_sites(sequence, binding_site='AGGAGG'):
    """
    Adds specific miRNA binding sites to a RNA sequence.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be modified.
    binding_site : str, optional
        The binding site to be added to the sequence (default is 'AGGAGG').

    Returns
    -------
    str
        The modified RNA sequence with the binding sites added.
    """
    binding_site_length = len(binding_site)
    return re.sub(f'(?=([A-Z]{{{binding_site_length}}}))', f'{binding_site}', sequence)


def convert_sequence_format(sequence, input_format='string', output_format='FASTA'):
    """
    Converts a RNA sequence from one format to another.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be converted.
    input_format : str, optional
        The format of the input sequence, either 'string' or 'FASTA' (default is 'string').
    output_format : str, optional
        The desired format of the output sequence, either 'string' or 'FASTA' (default is 'FASTA').

    Returns
    -------
    str
        The converted RNA sequence in the desired output format.
    """
    if input_format not in ['string', 'FASTA']:
        raise ValueError("Invalid input format. Choose either 'string' or 'FASTA'.")
    if output_format not in ['string', 'FASTA']:
        raise ValueError("Invalid output format. Choose either 'string' or 'FASTA'.")

    if input_format == 'string' and output_format == 'FASTA':
        return f'>sequence\n{sequence}'
    elif input_format == 'FASTA' and output_format == 'string':
        return sequence.split('\n')[1]
    else:
        return sequence


def export_sequence(sequence, output_file, output_format='FASTA'):
    """
    Exports a RNA sequence to a file in a specified format.

    Parameters
    ----------
    sequence : str
        The RNA sequence to a file in a specified format.

    Parameters
    ----------
    sequence : str
        The RNA sequence to be exported.
    output_file : str
        The path of the file to which the sequence will be exported.
    output_format : str, optional
        The format of the exported sequence, either 'string' or 'FASTA' (default is 'FASTA').

    Returns
    -------
    None
    """
    if output_format not in ['string', 'FASTA']:
        raise ValueError("Invalid output format. Choose either 'string' or 'FASTA'.")
    with open(output_file, 'w') as f:
        f.write(sequence)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A tool for analyzing and modifying RNA sequences.')
    parser.add_argument('sequence', type=str, help='The RNA sequence to be processed.')
    parser.add_argument('-a', '--analysis', action='store_true', help='Perform analysis on the sequence.')
    parser.add_argument('-m', '--modification', action='store_true', help='Perform modification on the sequence.')
    parser.add_argument('-e', '--export', type=str, help='Export the processed sequence to a file.')
    args = parser.parse_args()

    sequence = args.sequence
    if args.analysis:
        gc = gc_content(sequence)
        codons = codon_usage(sequence)
        splice_sites = splice_site_prediction(sequence)
        print(f'GC content: {gc}')
        print(f'Codon usage: {codons}')
        print(f'Splice sites: {splice_sites}')
    if args.modification:
        no_introns = remove_introns(sequence)
        binding_sites = add_mirna_binding_sites(sequence)
        fasta_format = convert_sequence_format(sequence, output_format='FASTA')
        print(f'Sequence without introns: {no_introns}')
        print(f'Sequence with binding sites: {binding_sites}')
        print(f'Sequence in FASTA format: {fasta_format}')
    if args.export:
        export_sequence(sequence, args.export)




