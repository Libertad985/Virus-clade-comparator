import os
from Bio import Align
from Bio import SeqIO


def pairwise_compare(reference, sample):
    """Compare two sequences and identify differences."""
    aligner = Align.PairwiseAligner(scoring="blastn")
    alignments = aligner.align(reference.seq, sample.seq)
    differences = set()
    for pos in range(len(alignments[0][1])):
        reference_base = alignments[0][:, pos][0]
        sample_base = alignments[0][:, pos][1]
        if reference_base == '-':
            differences.add((pos, 'ins', sample_base))
        elif sample_base == '-':
            differences.add((pos, 'del', reference_base))
        elif reference_base != sample_base:
            differences.add((pos, reference_base, sample_base))
    return differences


def main():
    # Set paths
    path_to_workdir = "<path_to_workdir>"
    reference_path = os.path.join(path_to_workdir, "<path_to_reference_file_in_fasta_format>")
    clade_path = os.path.join(path_to_workdir, "<path_to_clade_multifasta>")
    
    # Read reference sequence
    reference = SeqIO.read(reference_path, "fasta")
    
    # Parse clade sequences
    clade_sequences = SeqIO.parse(clade_path, "fasta")
    
    # Compare sequences
    intersection = None
    for sequence in clade_sequences:
        sequence.seq = str(sequence.seq).upper()
        differences = pairwise_compare(reference, sequence)
        if intersection is None:
            intersection = differences
        else:
            intersection = intersection.intersection(differences)
    
    # Print differences
    print("Differences:")
    for difference in intersection:
        print(difference)
if __name__ == "__main__":
    main()