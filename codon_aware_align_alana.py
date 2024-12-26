import argparse
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

def read_input_seqs(in_fasta):
    seq_records = SeqIO.parse(in_fasta, 'fasta')
    nucl_seqs = {}
    for rec in seq_records:
        nucl_seqs[rec.id] = rec
    return nucl_seqs

def write_aa_file(nucl_seqs, output_file):
    records = [rec.translate(table='Bacterial', id=rec.id, description=rec.description) for rec in nucl_seqs.values()]
    SeqIO.write(records, output_file, 'fasta')

def back_translate(aa_aln, original_seqs):
    nucl_aln = []
    for rec in aa_aln:
        num_gaps = 0
        aln_seq = []
        nucl_seq = original_seqs[rec.id]
        for i in range(len(rec.seq)):
            if rec[i] == '-':
                aln_seq.append('---')
                num_gaps += 1
            else:
                idx = 3 * (i - num_gaps)
                aln_seq.append(str(nucl_seq[idx:(idx + 3)].seq))
        nucl_aln.append(SeqRecord(Seq(''.join(aln_seq), generic_dna), id=rec.id, description=rec.description))
    return MultipleSeqAlignment(nucl_aln)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_fasta', type=str, help='Path to input FASTA file (must have .fasta extension).')
    parser.add_argument('-m', '--muscle_bin', type=str, default='muscle', help='Path to MUSCLE binary program.')
    parser.add_argument('-o', '--output_file', type=str, default=None, help='Output file for nucleotide alignments.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print MUSCLE output.')
    args = parser.parse_args()

    if args.output_file is None:
        f_out = args.in_fasta.replace('.fasta', '_aln.fasta')
    else:
        f_out = args.output_file

    nucl_seqs = read_input_seqs(args.in_fasta)
    f_aa_seqs = args.in_fasta.replace('.fasta', '.faa')
    write_aa_file(nucl_seqs, f_aa_seqs)

    # Perform alignment
    f_aa_aln = f_out.replace('.fasta', '.faa')
    muscle_stdout = subprocess.run([args.muscle_bin, '-in', f_aa_seqs, '-out', f_aa_aln], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    if args.verbose:
        print(muscle_stdout)

    aa_aln = AlignIO.read(f_aa_aln, 'fasta')
    nucl_aln = back_translate(aa_aln, nucl_seqs)
    AlignIO.write(nucl_aln, f_out, 'fasta')

