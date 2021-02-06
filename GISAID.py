import sys
import parasail
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import csv
import pandas

# system waits for user to import a file
# input_file = sys.stdin
input_file = "C:\Research\Spike\GISAID\substitution_and_indel_mutation_finder\spikeprot0129.fasta" # UNCOMMENT FOR OFFICIAL 
# input_file = "smalltest.fasta"

# Define the reference sequence: Below is the actual SARS-CoV-2 Protein obtained here: https://www.ncbi.nlm.nih.gov/protein/YP_009724390 
refseq_record = 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'
# this is just a simple sequence for testing
# refseq_record = 'ithinkpuppiesareverycoolHHHHHHHHHHHHHHHHHHHIIIIVWEOGBWNEO' # all of the sequences we want to compare to. 

# all of the sequences we want to compare to. 
seq_records = SeqIO.parse(input_file, 'fasta')

GAP_OPEN_PENALTY = 5
GAP_EXTEND_PENALTY = 2

def align(sequence1: str, sequence2: str): 
    """ Run the local pairwise alignment of two strings and return alignment data. """
    result = parasail.nw_trace(sequence1, 
                               sequence2, 
                               GAP_OPEN_PENALTY, 
                               GAP_EXTEND_PENALTY, 
                               parasail.blosum80)
    cigar_text = result.cigar.decode
    result = result.traceback.ref, result.traceback.query, result.score
    return result

def write_to_csv(mutations):
    df = pandas.DataFrame(mutations, columns=["type","mutations"])
    df.to_csv('GISAID.csv', chunksize=500000)

def main():
    count = 0
    insertions = []
    deletions = []
    substitutions = []

    for seq_record in seq_records:
        if len(seq_record.seq) < 800:
            None
        else:
            seq_aligned = list(align(str(seq_record.seq).upper(), refseq_record.upper()))
            # print(seq_aligned[0], seq_aligned[2], seq_aligned[1], sep='\n')
            aligned = [seq_aligned[0], seq_aligned[1]] 
            count += 1
            if count % 1000 == 0:
                print (count)
            for i in range(len(aligned[0])):
                nt1 = aligned[0]
                nt1 = nt1[i]
                nt2 = aligned[1]
                nt2 = nt2[i]
                
                if nt1 != nt2:
                    if nt2 == "-":
                        deletion = ("deletion", f"{nt1}{i+1}-")
                        deletions.append(deletion)
                    elif nt1 == "-":
                        insertion = ("insertion", f"{nt2}{i+1}{nt2}")
                        insertions.append(insertion)
                    elif nt1 == "*" or nt2 == "*":
                        None
                    else:
                        substitution = ("substitution", f"{nt1}{i+1}{nt2}")
                        substitutions.append(substitution)
                

    mutations = substitutions + insertions + deletions
    print(len(mutations))
    
    write_to_csv(mutations)
    
    
main()

