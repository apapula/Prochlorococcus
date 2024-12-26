import argparse
import numpy as np
import pandas as pd
import pickle
from Bio import AlignIO
from Bio.Seq import Seq


class CodonSitesDict:
    def __init__(self):
        self.initialize_synonymous_fraction_dict()

    def initialize_synonymous_fraction_dict(self):
        '''
        Manually define fraction of synonymous substitutions at each
        codon site from bacterial genetic code (table 11).
        '''
        self.codon_fraction_synonymous = {}
        self.codon_fraction_synonymous['TTT'] = np.array([0, 0, 1./3]) # Phenylalanine (F)
        self.codon_fraction_synonymous['TTC'] = np.array([0, 0, 1./3]) # Phenylalanine (F)
        self.codon_fraction_synonymous['TTA'] = np.array([1./3, 0, 1./3]) # Leucine (L)
        self.codon_fraction_synonymous['TTG'] = np.array([1./3, 0, 1./3]) # Leucine (L)
        self.codon_fraction_synonymous['TCT'] = np.array([0, 0, 1]) # Serine (S)
        self.codon_fraction_synonymous['TCC'] = np.array([0, 0, 1]) # Serine (S)
        self.codon_fraction_synonymous['TCA'] = np.array([0, 0, 1]) # Serine (S)
        self.codon_fraction_synonymous['TCG'] = np.array([0, 0, 1]) # Serine (S)
        self.codon_fraction_synonymous['TAT'] = np.array([0, 0, 1./3]) # Tyrosine (Y)
        self.codon_fraction_synonymous['TAC'] = np.array([0, 0, 1./3]) # Tyrosine (Y)
        self.codon_fraction_synonymous['TAA'] = np.array([0, 1./3, 1./3]) # Stop (Ochre)
        self.codon_fraction_synonymous['TAG'] = np.array([0, 0, 1./3]) # Stop (Amber)
        self.codon_fraction_synonymous['TGT'] = np.array([0, 0, 1./3]) # Cysteine (C)
        self.codon_fraction_synonymous['TGC'] = np.array([0, 0, 1./3]) # Cysteine (C)
        self.codon_fraction_synonymous['TGA'] = np.array([0, 1./3, 0]) # Stop (Opal)
        self.codon_fraction_synonymous['TGG'] = np.array([0, 0, 0]) # Tryptophan (W)

        self.codon_fraction_synonymous['CTT'] = np.array([0, 0, 1]) # Leucine (L)
        self.codon_fraction_synonymous['CTC'] = np.array([0, 0, 1]) # Leucine (L)
        self.codon_fraction_synonymous['CTA'] = np.array([1./3, 0, 1]) # Leucine (L)
        self.codon_fraction_synonymous['CTG'] = np.array([1./3, 0, 1]) # Leucine (L)
        self.codon_fraction_synonymous['CCT'] = np.array([0, 0, 1]) # Proline (P)
        self.codon_fraction_synonymous['CCC'] = np.array([0, 0, 1]) # Proline (P)
        self.codon_fraction_synonymous['CCA'] = np.array([0, 0, 1]) # Proline (P)
        self.codon_fraction_synonymous['CCG'] = np.array([0, 0, 1]) # Proline (P)
        self.codon_fraction_synonymous['CAT'] = np.array([0, 0, 1./3]) # Histidine (H)
        self.codon_fraction_synonymous['CAC'] = np.array([0, 0, 1./3]) # Histidine (H)
        self.codon_fraction_synonymous['CAA'] = np.array([0, 0, 1./3]) # Glutamine (Q)
        self.codon_fraction_synonymous['CAG'] = np.array([0, 0, 1./3]) # Glutamine (Q)
        self.codon_fraction_synonymous['CGT'] = np.array([0, 0, 1]) # Arginine (R)
        self.codon_fraction_synonymous['CGC'] = np.array([0, 0, 1]) # Arginine (R)
        self.codon_fraction_synonymous['CGA'] = np.array([1./3, 0, 1]) # Arginine (R)
        self.codon_fraction_synonymous['CGG'] = np.array([1./3, 0, 1]) # Arginine (R)

        self.codon_fraction_synonymous['ATT'] = np.array([0, 0, 2./3]) # Isoleucine (I)
        self.codon_fraction_synonymous['ATC'] = np.array([0, 0, 2./3]) # Isoleucine (I)
        self.codon_fraction_synonymous['ATA'] = np.array([0, 0, 2./3]) # Isoleucine (I)
        self.codon_fraction_synonymous['ATG'] = np.array([0, 0, 0]) # Methionine (M)
        self.codon_fraction_synonymous['ACT'] = np.array([0, 0, 1]) # Threonine (T)
        self.codon_fraction_synonymous['ACC'] = np.array([0, 0, 1]) # Threonine (T)
        self.codon_fraction_synonymous['ACA'] = np.array([0, 0, 1]) # Threonine (T)
        self.codon_fraction_synonymous['ACG'] = np.array([0, 0, 1]) # Threonine (T)
        self.codon_fraction_synonymous['AAT'] = np.array([0, 0, 1./3]) # Asparagine (N)
        self.codon_fraction_synonymous['AAC'] = np.array([0, 0, 1./3]) # Asparagine (N)
        self.codon_fraction_synonymous['AAA'] = np.array([0, 0, 1./3]) # Lysine (K)
        self.codon_fraction_synonymous['AAG'] = np.array([0, 0, 1./3]) # Lysine (K)
        self.codon_fraction_synonymous['AGT'] = np.array([0, 0, 1./3]) # Serine (S)
        self.codon_fraction_synonymous['AGC'] = np.array([0, 0, 1./3]) # Serine (S)
        self.codon_fraction_synonymous['AGA'] = np.array([1./3, 0, 1./3]) # Arginine (R)
        self.codon_fraction_synonymous['AGG'] = np.array([1./3, 0, 1./3]) # Arginine (R)

        self.codon_fraction_synonymous['GTT'] = np.array([0, 0, 1]) # Valine (V)
        self.codon_fraction_synonymous['GTC'] = np.array([0, 0, 1]) # Valine (V)
        self.codon_fraction_synonymous['GTA'] = np.array([0, 0, 1]) # Valine (V)
        self.codon_fraction_synonymous['GTG'] = np.array([0, 0, 1]) # Valine (V)
        self.codon_fraction_synonymous['GCT'] = np.array([0, 0, 1]) # Alanine (A)
        self.codon_fraction_synonymous['GCC'] = np.array([0, 0, 1]) # Alanine (A)
        self.codon_fraction_synonymous['GCA'] = np.array([0, 0, 1]) # Alanine (A)
        self.codon_fraction_synonymous['GCG'] = np.array([0, 0, 1]) # Alanine (A)
        self.codon_fraction_synonymous['GAT'] = np.array([0, 0, 1./3]) # Aspartic acid (D)
        self.codon_fraction_synonymous['GAC'] = np.array([0, 0, 1./3]) # Aspartic acid (D)
        self.codon_fraction_synonymous['GAA'] = np.array([0, 0, 1./3]) # Glutamic acid (E)
        self.codon_fraction_synonymous['GAG'] = np.array([0, 0, 1./3]) # Glutamic acid (E)
        self.codon_fraction_synonymous['GGT'] = np.array([0, 0, 1]) # Glycine (G)
        self.codon_fraction_synonymous['GGC'] = np.array([0, 0, 1]) # Glycine (G)
        self.codon_fraction_synonymous['GGA'] = np.array([0, 0, 1]) # Glycine (G)
        self.codon_fraction_synonymous['GGG'] = np.array([0, 0, 1]) # Glycine (G)

    def get_synonymous_substitutions_fraction(self, codon):
        return self.codon_fraction_synonymous[codon]

    def get_site_degeneracies(self, codon):
        f_synonymous = self.get_synonymous_substitutions_fraction(codon)
        degeneracies = []
        for f in f_synonymous:
            if f == 0:
                degeneracies.append('1D')
            elif f > 0 and f < 2./3:
                degeneracies.append('2D')
            elif f > 1./3 and f < 1:
                degeneracies.append('3D')
            else:
                degeneracies.append('4D')
        return degeneracies

    def is_codon(self, string):
        return string in self.codon_fraction_synonymous


def calculate_pairwise_pNpS(aln):
    sag_ids = [record.id for record in aln]
    pN_df = pd.DataFrame(index=sag_ids, columns=sag_ids)
    pS_df = pd.DataFrame(index=sag_ids, columns=sag_ids)
    for i in range(len(aln)):
        seq1 = aln[i]
        pN_df.at[seq1.id, seq1.id] = 0
        pS_df.at[seq1.id, seq1.id] = 0
        for j in range(i):
            seq2 = aln[j]

            pN, pS = calculate_ng86_pdist(seq1, seq2)
            pN_df.at[seq1.id, seq2.id] = pN
            pN_df.at[seq2.id, seq1.id] = pN
            pS_df.at[seq1.id, seq2.id] = pS
            pS_df.at[seq2.id, seq1.id] = pS
    return pN_df, pS_df


def calculate_ng86_pdist(seq1, seq2):
    '''
    Uses method II form Nei & Gojobori (1986) to calculate fraction of substituions
    between two sequences.
    '''

    codon_dict = CodonSitesDict()

    num_codons = 0
    num_nongapped_codons = 0
    S = 0 # num synonymous sites
    N = 0 # num non-synonymous sites
    m_s = np.zeros(3) # num substitution sites at codons with single difference
    m_ss = np.zeros(3) # num synonymous sites at codons with single difference
    m_m = np.zeros(3) # num substitution sites at codons with multiple differences

    if len(seq1) == len(seq2):
        # Read sequence codon by codon
        for i in range(0, len(seq1), 3):
            codon1 = seq1.seq[i:i + 3]
            codon2 = seq2.seq[i:i + 3]

            if '-' in codon1 or '-' in codon2:
                #print(i, codon1, codon2)
                # Ignore codon if gaps are present
                continue
            elif 'N' in codon1 or 'N' in codon2:
                # Ignore codons with 'N's
                continue
            elif len(codon1) == 3 and len(codon2):
                # Gaps that are not multiples of 3 can lead to bad reading frame
                # Ignoring for now. TODO: find workaround
                num_nongapped_codons += 1
                f1 = codon_dict.get_synonymous_substitutions_fraction(codon1)
                f2 = codon_dict.get_synonymous_substitutions_fraction(codon2)
                site_diffs = np.array(codon1) != np.array(codon2)
                num_diffs = np.sum(site_diffs, dtype=np.int64)
                if num_diffs == 1:
                    m_s[site_diffs] += 1
                    m_ss[site_diffs] += (f1[site_diffs] + f2[site_diffs]) / 2
                    #print(codon1, codon2)
                elif num_diffs > 1:
                    m_m[site_diffs] += 1
                    #print(site_diffs, codon1, codon2)
                S += np.sum((f1 + f2) / 2)
                N += 3 - np.sum((f1 + f2) / 2)

    #print(m_ss, m_s, num_nongapped_codons)
    pi_s = m_ss / (m_s + (m_s == 0)) # avoid division by zero
    S_d = np.dot(m_s + m_m, pi_s)
    N_d = np.dot(m_s + m_m, 1 - pi_s)

    return N_d / N, S_d / S


def read_mafft_alignment(f_in, file_format='fasta', alphabet='generic_dna'):
    aln = AlignIO.read(f_in, file_format)
    for rec in aln:
        rec.seq = rec.seq.upper()
    return aln


def read_alignment(f_in, file_format='fasta_unambiguous', seq_type='nucl'):
    '''
    Assumes input is FASTA. Converts alphabet to upper case if necessary and replaces
        all 'N' with '-'.
    '''

    aln = AlignIO.read(f_in, 'fasta')
    if file_format == 'mafft':
        for record in aln:
            record.seq = record.seq.upper()

    if file_format != 'fasta_unambiguous':
        for record in aln:
            record.seq = replace_ambiguous_chars(record.seq, seq_type)

    return aln


def replace_ambiguous_chars(in_seq, seq_type):
    if seq_type == 'nucl':
        alphabet = ['A', 'T', 'G', 'C']
    else:
        alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    seq_letters = set(str(in_seq))
    if seq_letters - set(alphabet):
        # Replace ambiguous letters with gaps
        seq_arr = np.array(in_seq)
        replacement_idx = [n not in alphabet for n in seq_arr]
        seq_arr[replacement_idx] = '-'
        out_seq = Seq(''.join(seq_arr))
    else:
        out_seq = in_seq
    return out_seq



def test_pNpS():
    f_test = '../results/tests/YSG_0947_alpha_block1_seqs_v1.fna'
    aln = read_alignment(f_test)
    print(aln)

    pN_segment, pS_segment = calculate_pairwise_pNpS(aln)
    print(pN_segment)
    print(pS_segment)

    f_test = '../results/tests/YSG_0947_alpha_block1_seqs_v2.fna'
    aln = read_alignment(f_test)
    pN_segment, pS_segment = calculate_pairwise_pNpS(aln)
    print(aln)
    print(pN_segment)
    print(pS_segment)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Alignment file in FASTA format.')
    parser.add_argument('-f', '--file_format', default='default', help='FASTA file format ["default", "mafft"].')
    parser.add_argument('-o', '--output_file', default=None)
    parser.add_argument('-t', '--test', action='store_true')
    parser.add_argument('--output_format', default='pickle', help='Output file format ["pickle", "tsv"]. If "tsv", output_file should end in ".tsv"')
    args = parser.parse_args()

    if args.test == True:
        test_pNpS()
    else:
        if args.file_format == 'default':
            aln = read_alignment(args.input_file)
        else:
            aln = read_mafft_alignment(args.input_file)
        pN, pS = calculate_pairwise_pNpS(aln)

        if args.output_file is not None:
            if args.file_format == 'pickle':
                with open(args.output_file, 'wb') as fout:
                    pickle.dump({'pN':pN, 'pS':pS}, fout)
            elif args.output_format == 'tsv':
                fN = args.output_file.replace('.tsv', '_pN.tsv')
                pN.to_csv(fN, sep='\t')
                fS = args.output_file.replace('.tsv', '_pS.tsv')
                pS.to_csv(fS, sep='\t')
        else:
            print(f'pN:\n{pN}\n\n')
            print(f'pS:\n{pS}')

