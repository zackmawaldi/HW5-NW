# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    non_hs_seqs = [gg_seq, mm_seq, br_seq, tt_seq]
    aligned_seqs = []

    for seq in non_hs_seqs:
        nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
        nw.align(hs_seq, seq)
        aligned_seqs.append( (nw.alignment_score, nw.seqA_align, nw.seqB_align) )
    
    aligned_seqs = sorted( aligned_seqs , reverse=True)

    for alignment_tuple in aligned_seqs:
        print(alignment_tuple[1])
        print(alignment_tuple[2])
        print('Score:', alignment_tuple[0])
        print()
    

if __name__ == "__main__":
    main()
