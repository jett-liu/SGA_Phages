#!/home/linking/.pyenv/versions/3.8.2/bin/python
# Author: Lin-Xing Chen, UC Berkeley

# CRISPR.py
# Version 0.4
# 7.30.2021


import re
import os
import argparse
import time
from time import strftime
import pysam
from Bio import SeqIO
from Bio.Seq import reverse_complement
from collections import defaultdict


parser = argparse.ArgumentParser(description="This script gets all spacers from both paired reads when providing assembled sequences and the read mapping file.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-f", "--fasta", type=str, help="The scaffolds file", required=True)
requiredNamed.add_argument("-m", "--mapping", type=str, help="Scaffolds mapping sam/bam file.", required=True)
parser.add_argument("-d", "--database", type=str, help="cas protein HMM database provided by user, otherwise the TIGRFAM Cas HMM will be used")
parser.add_argument("-c", "--cpu", type=int, default=16, help="Number of cpu for hmmscan (default=16).")
parser.add_argument("-n", "--process", type=int, default=16, help="Number of processes for blastn (default=16).")
parser.add_argument("-o", "--output", type=str, help="The output folder.")
args = parser.parse_args()



# TIGRFAMs IDs of cas proteins
tigrID2cas = {'TIGR00287': 'cas1', 'TIGR00372': 'cas4', 'TIGR01573': 'cas2', 'TIGR01587': 'cas3_core',
               'TIGR01595': 'cas_CT1132', 'TIGR01596': 'cas3_HD', 'TIGR01863': 'cas_Csd1', 'TIGR01865': 'cas_Csn1',
               'TIGR01866': 'cas_Csn2', 'TIGR01868': 'casD_Cas5e', 'TIGR01869': 'casC_Cse4',
               'TIGR01870': 'cas_TM1810_Csm2', 'TIGR01873': 'cas_CT1978', 'TIGR01874': 'cas_cas5a',
               'TIGR01875': 'cas_MJ0381', 'TIGR01876': 'cas_Cas5d', 'TIGR01877': 'cas_cas6', 'TIGR01878': 'cas_Csa5',
               'TIGR01881': 'cas_Cmr5', 'TIGR01884': 'cas_HTH', 'TIGR01888': 'cas_cmr3', 'TIGR01894': 'cas_TM1795_cmr1',
               'TIGR01895': 'cas_Cas5t', 'TIGR01896': 'cas_AF1879', 'TIGR01897': 'cas_MJ1666',
               'TIGR01898': 'cas_TM1791_cmr6', 'TIGR01899': 'cas_TM1807_csm5', 'TIGR01903': 'cas5_csm4',
               'TIGR01907': 'casE_Cse3', 'TIGR01908': 'cas_CXXC_CXXC', 'TIGR01914': 'cas_Csa4',
               'TIGR02165': 'cas5_6_GSU0054', 'TIGR02221': 'cas_TM1812', 'TIGR02547': 'casA_cse1',
               'TIGR02548': 'casB_cse2', 'TIGR02549': 'CRISPR_DxTHG', 'TIGR02556': 'cas_TM1802',
               'TIGR02562': 'cas3_yersinia', 'TIGR02563': 'cas_Csy4', 'TIGR02564': 'cas_Csy1',
               'TIGR02565': 'cas_Csy2', 'TIGR02566': 'cas_Csy3', 'TIGR02570': 'cas7_GSU0053',
               'TIGR02577': 'cas_TM1794_Cmr2', 'TIGR02578': 'cas_TM1811_Csm1', 'TIGR02579': 'cas_csx3',
               'TIGR02580': 'cas_RAMP_Cmr4', 'TIGR02581': 'cas_cyan_RAMP', 'TIGR02582': 'cas7_TM1809',
               'TIGR02583': 'DevR_archaea', 'TIGR02584': 'cas_NE0113', 'TIGR02585': 'cas_Cst2_DevR',
               'TIGR02586': 'cas5_cmx5_devS', 'TIGR02589': 'cas_Csd2', 'TIGR02590': 'cas_Csh2',
               'TIGR02591': 'cas_Csh1', 'TIGR02592': 'cas_Cas5h', 'TIGR02593': 'CRISPR_cas5', 'TIGR02619': 'TIGR02619',
               'TIGR02620': 'cas_VVA1548', 'TIGR02621': 'cas3_GSU0051', 'TIGR02670': 'cas_csx8',
               'TIGR02671': 'cas_csx9', 'TIGR02672': 'cas_csm6', 'TIGR02674': 'cas_cyan_RAMP_2',
               'TIGR02682': 'cas_csx11', 'TIGR02710': 'TIGR02710', 'TIGR02807': 'cas6_cmx6', 'TIGR03031': 'cas_csx12',
               'TIGR03114': 'cas8u_csf1', 'TIGR03115': 'cas7_csf2', 'TIGR03116': 'cas5_csf3', 'TIGR03117': 'cas_csf4',
               'TIGR03157': 'cas_Csc2', 'TIGR03158': 'cas3_cyano', 'TIGR03159': 'cas_Csc1', 'TIGR03174': 'cas_Csc3',
               'TIGR03485': 'cas_csx13_N', 'TIGR03486': 'cas_csx13_C', 'TIGR03487': 'cas_csp2',
               'TIGR03488': 'cas_Cas5p', 'TIGR03489': 'cas_csp1', 'TIGR03637': 'cas1_YPEST', 'TIGR03638': 'cas1_ECOLI',
               'TIGR03639': 'cas1_NMENI', 'TIGR03640': 'cas1_DVULG', 'TIGR03641': 'cas1_HMARI',
               'TIGR03642': 'cas_csx14', 'TIGR03876': 'cas_csaX', 'TIGR03983': 'cas1_MYXAN', 'TIGR03984': 'TIGR03984',
               'TIGR03985': 'TIGR03985', 'TIGR03986': 'TIGR03986', 'TIGR04093': 'cas1_CYANO',
               'TIGR04106': 'cas8c_GSU0052', 'TIGR04113': 'cas_csx17', 'TIGR04328': 'cas4_PREFRAN',
               'TIGR04329': 'cas1_PREFRAN', 'TIGR04330': 'cas_Cpf1', 'TIGR04413': 'MYXAN_cmx8',
               'TIGR04423': 'casT3_TIGR04423', 'COG3513': 'COG3513_Cas9', 'cd09643': 'cd09643_Cas9', 'cd09704': 'cd09704_Cas9',
               'icity0002': 'icity0002_Cas9', 'mkCas0193': 'mkCas0193_Cas9'}


#
def get_gene_id(scaffold_id, align_start, align_end):
    for gene in gene2loc.keys():
        if gene.rsplit('_', 1)[0] == scaffold_id and int(gene2loc[gene][0]) <= int(align_start) and int(gene2loc[gene][1]) >= int(align_end):
            return gene
        else:
            continue
    return "Position: " + align_start + '-' + align_end


#
def is_nucleotide(item):
    total = 0
    bases = ['A', 'T', 'C', 'G']
    for base in item:
        if base not in bases:
            total += 1
        else:
            pass
    return total == 0


#
def get_the_longest(List):
    longest = 0
    longest_seq = ''
    for item in List:
        header = item.rsplit('_region_', 1)[0]
        if len(header2seq[header]) >= longest and item in cas_sequences:
            longest = len(header2seq[header])
            longest_seq = item
        else:
            pass
    return longest_seq


#
def one_has_cas(direct_repeat):
    total = 0
    target = ''
    for item in dr2region[direct_repeat]:
        if item in cas_sequences:
            total += 1
            target = item
    if total == 1:
        return target
    elif total > 1:
        return get_the_longest(dr2region[direct_repeat])
    else:
        return "nothing"


# the current path
pwd = os.getcwd()
fasta_file = os.path.abspath(args.fasta)
mapping_file = os.path.abspath(args.mapping)

# Check if output directory exists or not
fasta_name = args.fasta.strip().split('/')[-1]
if args.output:
    if os.path.exists('{0}'.format(args.output)):
        print('Warning: Output directory exists, please rename the existed one!')
        exit()
    else:
        os.mkdir('{0}'.format(args.output))
        os.chdir('{0}'.format(args.output))
else:
    if os.path.exists("{0}.crispr".format(fasta_name)):
        print('Warning: Output directory exists, please rename the existed one!')
        exit()
    else:
        os.mkdir("{0}.crispr".format(fasta_name))
        os.chdir("{0}.crispr".format(fasta_name))


# log file
log = open("log", 'w+')
star_num = int((130-len(' CRISPR analyses for {0} '.format(fasta_name)))/2)
print('*' * star_num + ' CRISPR analyses for {0} '.format(fasta_name) + '*' * star_num, end='\n', file=log, flush=True)


# Open scaffolds file
print('[{0}] [1/12] Reading sequence file ... '.format(strftime("%c")), end='', file=log, flush=True)
header2seq = {}
file_size = 0
with open('{0}'.format(fasta_file), 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
        header = str(record.id).strip()
        seq = str(record.seq)
        header2seq[header] = seq
        file_size += len(seq)
print('A total of {0} sequences were imported ...'.format(str(len(header2seq.keys()))), end='\n', file=log, flush=True)


# predict direct repeat region from scaffolds/contigs
print('[{0}] [2/12] Searching sequences for direct repeat region by Pilercr ...'.format(strftime("%c")), end='\n', file=log, flush=True)
if file_size < 500000000:
    os.system("/groups/banfield/users/linking/softwares/pilercr1.06/./pilercr -in {0} -out {1}.pilercrl -noinfo -seq {1}.pilercrl.seq".format(fasta_file, fasta_name))
else:  # if the fasta file is larger than 500 Mbp, it should be split into smaller ones, otherwise pilercr will not work
    accumulated_length = 0
    header2file = {}
    for header in header2seq.keys():
        accumulated_length += len(header2seq[header])
        header2file[header] = str(int(accumulated_length / 500000000))

    os.system('mkdir split_fasta_files')
    split_fasta = {}
    file_lists = []
    for header in header2file.keys():
        id = header2file[header]
        if id not in split_fasta:
            split_fasta[id] = open('split_fasta_files/{0}.{1}'.format(fasta_name, id), 'w')
            file_lists.append('{0}.{1}'.format(fasta_name, id))
            split_fasta[id].write('>' + header + '\n')
            split_fasta[id].write(header2seq[header] + '\n')
        else:
            split_fasta[id].write('>' + header + '\n')
            split_fasta[id].write(header2seq[header] + '\n')
    for id in split_fasta.keys():
        split_fasta[id].close()

    for file in file_lists:
        os.system("/groups/banfield/users/linking/softwares/pilercr1.06/./pilercr -in split_fasta_files/{0} -out split_fasta_files/{0}.pilercrl -noinfo "
                  "-seq split_fasta_files/{0}.pilercrl.seq".format(file))

    os.system('cat split_fasta_files/*pilercrl >{0}.pilercrl'.format(fasta_name))
    os.system('cat split_fasta_files/*pilercrl.seq >{0}.pilercrl.seq'.format(fasta_name))


# parse the Pilercr results
print('[{0}] [3/12] Paring Pilercr for direct repeat region ...'.format(strftime("%c")), end='\n', file=log, flush=True)
#
sequences_for_prodigal = set()
keys = []
scaff2array2start_repeat = {}
pilercrl_seq = open('{0}.pilercrl.seq'.format(fasta_name), 'r')
for line in pilercrl_seq.readlines():
    if line.startswith('>'):
        line = line.strip().rsplit('[',1)
        scaffold = line[0].split(' ')[0][1:]
        sequences_for_prodigal.add(scaffold)
        array = line[1].split(';')[0][5:]
        start = line[1].rsplit('=',1)[1][:-1]
        keys.append(scaffold)
        keys.append(array)
        if scaffold not in scaff2array2start_repeat.keys():
            scaff2array2start_repeat[scaffold] = {}
            scaff2array2start_repeat[scaffold][array] = [start]
        else:
            scaff2array2start_repeat[scaffold][array] = [start]
    else:
        pass
pilercrl_seq.close()

#
pilercrl = open('{0}.pilercrl'.format(fasta_name), 'r')
scaffold_lists = []
for line in pilercrl.readlines():
    while '  ' in line:
        line = line.replace('  ', ' ')
    line = line.strip().split(' ')
    if line[0].startswith('>'):
        scaffold_lists.append(line[0].strip()[1:])
    elif len(line) >= 8 and is_nucleotide(line[-1]) and line[-2] != '-' and line[-2] != '+':
        scaff2array2start_repeat[scaffold_lists[-1]][line[0]].append(str(int(line[2]) + int(line[3]) - 1))
        scaff2array2start_repeat[scaffold_lists[-1]][line[0]].append(line[-1])
    else:
        pass
pilercrl.close()


dr_region2repeat = {}
dr2region = {}

for scaffold in scaff2array2start_repeat.keys():
    for array in scaff2array2start_repeat[scaffold].keys():
        start = scaff2array2start_repeat[scaffold][array][0]
        end = scaff2array2start_repeat[scaffold][array][1]
        dr = scaff2array2start_repeat[scaffold][array][2]
        region_name = scaffold + '_region_' + start + '_' + end
        dr_region2repeat[region_name] = dr
        if dr not in dr2region.keys():
            if reverse_complement(dr) in dr2region.keys():
                dr2region[reverse_complement(dr)].append(region_name)
            else:
                dr2region[dr] = [region_name]
        else:
            dr2region[dr].append(region_name)


# for debug
print(dr2region, flush=True)


print('[{0}] [4/12] Extracting and predicting genes for sequences with direct repeat region ...'.format(strftime("%c")), end='\n', file=log, flush=True)
# extract sequences with direct repeat region(s)
with open('sequences.with.direct.repeat.fasta', 'w') as s:
    for name in sequences_for_prodigal:
        s.write('>' + name + '\n')
        s.write(header2seq[name] + '\n')
s.close()

# predict genes using Prodigal
os.system('prodigal -i {0} -a {0}.genes.faa -d {0}.genes.fna -m -p meta'.format('sequences.with.direct.repeat.fasta'))

# extract genes within 10k distance
for_cas = open('sequences.with.direct.repeat.predicted.genes.within.10k.faa', 'w')
with open('sequences.with.direct.repeat.fasta.genes.faa', 'r') as faa:
    for record in SeqIO.parse(faa, "fasta"):
        header = str(record.id).strip()
        seq = str(record.seq)
        for item in dr_region2repeat.keys():
            if header.rsplit('_', 1)[0] == item.rsplit('_region_')[0]:
                dr_region_start = int(item.split('_')[-2])
                dr_region_end = int(item.split('_')[-1])
                if dr_region_start - 10000 < int(record.description.strip().split(' # ')[1]) < dr_region_start:
                    for_cas.write('>' + header + '_region_' + str(dr_region_start) + '_' + str(dr_region_end) + '\n')
                    for_cas.write(seq + '\n')
                elif dr_region_end < int(record.description.strip().split(' # ')[2]) < dr_region_end + 10000:
                    for_cas.write('>' + header + '_region_' + str(dr_region_start) + '_' + str(dr_region_end) + '\n')
                    for_cas.write(seq + '\n')
                else:
                    pass
            else:
                pass
for_cas.close()
faa.close()


# run HMM search for Cas proteins
print('[{0}] [5/12] Predicting Cas proteins using HMM databases ...'.format(strftime("%c")), end='\n', file=log, flush=True)
if args.database:
    os.system("hmmsearch --domtblout cas.proteins.hmm.txt --cpu {1} {2} {3}".
              format(fasta_name, args.cpu, args.database, 'sequences.with.direct.repeat.predicted.genes.within.10k.faa'))
else:
    os.system("hmmsearch --cut_nc --domtblout cas.proteins.hmm.txt --cpu {1} /groups/banfield/users/linking/databases/cas9.and.tigrfam.cas.hmm {2}".
              format(fasta_name, args.cpu, 'sequences.with.direct.repeat.predicted.genes.within.10k.faa'))


# parse the HMM search results using cath
os.system('/home/linking/softwares/./cath-resolve-hits.ubuntu14.04 --input-format hmmer_domtblout cas.proteins.hmm.txt >cas.proteins.hmm.cath.txt')

# parse the HMM search results
print('[{0}] [6/12] Parsing cas proteins HMM search results ...'.format(strftime("%c")), end='\n', file=log, flush=True)
cas_hmm_parsed = open("cas.proteins.hmm.cath.parsed.txt", 'w')
cas_hmm_parsed.write('Gene' + '\t' + 'TIGR_ID' + '\t' + 'Cas' + '\n')

cas_sequences = set()
dr_region_2_cas = {}

cas_hmm = open("cas.proteins.hmm.cath.txt", 'r')
for line in cas_hmm.readlines():
    if not line.startswith('#'):
        line = line.strip().split(' ')
        cas_hmm_parsed.write(line[0] + '\t' + line[1] + '\t' + tigrID2cas[line[1]] + '\n')
        sequence_name = line[0].rsplit('_region_', 1)[0].rsplit('_', 1)[0] + '_region_' + line[0].rsplit('_region_', 1)[1]
        if sequence_name not in dr_region_2_cas.keys():
            cas_sequences.add(sequence_name)
            dr_region_2_cas[sequence_name] = [tigrID2cas[line[1]]]
        else:
            cas_sequences.add(sequence_name)
            dr_region_2_cas[sequence_name].append(tigrID2cas[line[1]])
cas_hmm.close()
cas_hmm_parsed.close()


# link sequence with direct repeat but without Cas with sequence with both direct repeat and Cas via repeat
print('[{0}] [7/12] Checking sequences that sharing direct repeat ...'.format(strftime("%c")), end='\n', file=log, flush=True)
link_via_dr = {}
for dr in dr2region.keys():
    if one_has_cas(dr) != 'nothing':
        for region in dr2region[dr]:
            scaffold = region.rsplit('_region_', 1)[0]
            if scaffold not in link_via_dr.keys():
                link_via_dr[scaffold] = set()
                link_via_dr[scaffold].add(one_has_cas(dr))
            else:
                link_via_dr[scaffold].add(one_has_cas(dr))
    else:
        pass

for key in link_via_dr.keys():
    while None in link_via_dr[key]:
        link_via_dr[key].remove(None)  # None was included in the list due to the function of get_the_longest

# for debug
print(link_via_dr, flush=True)

# extract reads from *sam/*bam file
print('[{0}] [8/12] Extracting reads from read mapping file ...'.format(strftime("%c")), end='\n', file=log, flush=True)

dr_region_2_sequence = {}
for scaffold in link_via_dr.keys():
    for item in link_via_dr[scaffold]:
        #if item in dr_region2repeat.keys():
        if dr_region2repeat[item] in header2seq[scaffold]:
            if item not in dr_region_2_sequence.keys():
                dr_region_2_sequence[item] = [header2seq[scaffold]]
            else:
                dr_region_2_sequence[item].append(header2seq[scaffold])
        elif dr_region2repeat[item] in reverse_complement(header2seq[scaffold]):
            if item not in dr_region_2_sequence.keys():
                dr_region_2_sequence[item] = [reverse_complement(header2seq[scaffold])]
            else:
                dr_region_2_sequence[item].append(reverse_complement(header2seq[scaffold]))
        else:
            pass
        #else:
            #pass


dr_region_2_read = {}
map_file = pysam.AlignmentFile(mapping_file)
for line in map_file:
    if line.reference_name in link_via_dr.keys():
        for item in link_via_dr[line.reference_name]:
            if dr_region2repeat[item] in line.query_sequence:
                if item in dr_region_2_read.keys():
                    dr_region_2_read[item].append(line.query_sequence)
                else:
                    dr_region_2_read[item] = [line.query_sequence]
            elif dr_region2repeat[item] in reverse_complement(line.query_sequence):
                if item in dr_region_2_read.keys():
                    dr_region_2_read[item].append(reverse_complement(line.query_sequence))
                else:
                    dr_region_2_read[item] = [reverse_complement(line.query_sequence)]
            else:
                pass
    else:
        pass
map_file.close()

# if nothing is from the read
if dr_region_2_read == {}:
    print('Note: The sequence names in mapping file are different from those in your fasta file, no spacer was extracted from reads.', file=log, flush=True)
else:
    pass

# for debug
print(dr_region_2_sequence, flush=True)
print(dr_region_2_read, flush=True)


# extract spacers from both reads and sequences
print('[{0}] [9/12] Extracting spacers from reads and sequences ...'.format(strftime("%c")), end='\n', file=log, flush=True)
#
dr_region_2_half_repeat = {}
for dr_region in dr_region_2_sequence.keys():
    if len(dr_region2repeat[dr_region]) % 2 == 1:
        length = len(dr_region2repeat[dr_region]) + 1
        dr_region_2_half_repeat[dr_region] = [dr_region2repeat[dr_region][:int(length / 2)]]
        dr_region_2_half_repeat[dr_region].append(dr_region2repeat[dr_region][int(length / 2):])
    else:
        length = len(dr_region2repeat[dr_region])
        dr_region_2_half_repeat[dr_region] = [dr_region2repeat[dr_region][:int(length / 2)]]
        dr_region_2_half_repeat[dr_region].append(dr_region2repeat[dr_region][int(length / 2):])

# for debug
print(dr_region_2_half_repeat, flush=True)

#
dr_region_2_sequence_spacer = {}
for dr_region in dr_region_2_sequence.keys():
    for sequence in dr_region_2_sequence[dr_region]:
        if dr_region in dr_region_2_half_repeat.keys():
            for spacer in re.findall(dr_region_2_half_repeat[dr_region][1] + "(.+?)" + dr_region_2_half_repeat[dr_region][0], sequence):
                if 51 > len(spacer) > 15 and dr_region_2_half_repeat[dr_region][0] not in spacer and dr_region_2_half_repeat[dr_region][1] not in spacer:
                    if dr_region in dr_region_2_sequence_spacer.keys():
                        if spacer not in dr_region_2_sequence_spacer[dr_region]:
                            dr_region_2_sequence_spacer[dr_region].append(spacer)
                        else:
                            pass
                    else:
                        dr_region_2_sequence_spacer[dr_region]= [spacer]
                else:
                    pass
        else:
            pass


dr_region_2_read_spacer = {}
for dr_region in dr_region_2_read.keys():
    for read in dr_region_2_read[dr_region]:
        if dr_region in dr_region_2_half_repeat.keys():
            for spacer in re.findall(dr_region_2_half_repeat[dr_region][1] + "(.+?)" + dr_region_2_half_repeat[dr_region][0], read):
                if 51 > len(spacer) > 15 and dr_region_2_half_repeat[dr_region][0] not in spacer and dr_region_2_half_repeat[dr_region][1] not in spacer:
                    if dr_region in dr_region_2_read_spacer.keys():
                        if spacer not in dr_region_2_read_spacer[dr_region]:
                            dr_region_2_read_spacer[dr_region].append(spacer)
                        else:
                            pass
                    else:
                        dr_region_2_read_spacer[dr_region] = [spacer]
                else:
                    pass
        else:
            pass


# for debug
print(dr_region_2_sequence_spacer, flush=True)
print(dr_region_2_read_spacer, flush=True)


spacer_dic = {}
for dr_region in dr_region_2_sequence_spacer.keys():
    spacer_dic[dr_region] = []
    for spacer in dr_region_2_sequence_spacer[dr_region]:
        spacer_dic[dr_region].append(spacer)

    if dr_region in dr_region_2_read_spacer.keys():
        for spacer in dr_region_2_read_spacer[dr_region]:
            if spacer in spacer_dic[dr_region]:
                pass
            else:
                spacer_dic[dr_region].append(spacer)
    else:
        pass


spacer = open('spacers.fasta', 'w')
spacer_id_2_seq = {}
spacer_seq_2_dr_region = {}
for dr_region in spacer_dic.keys():
    for id in range(0, len(spacer_dic[dr_region])):
        if spacer_dic[dr_region][id] not in spacer_seq_2_dr_region.keys():
            spacer.write('>' + dr_region + '_spacer_' + str(id) + '\n')
            spacer.write(spacer_dic[dr_region][id] + '\n')
            spacer_id_2_seq[dr_region + '_spacer_' + str(id)] = spacer_dic[dr_region][id]
            spacer_seq_2_dr_region[spacer_dic[dr_region][id]] = dr_region
        else:
            if dr_region.rsplit('_region_', 1)[0] == spacer_seq_2_dr_region[spacer_dic[dr_region][id]].rsplit('_region_', 1)[0]:
                pass
            else:
                spacer.write('>' + dr_region + '_spacer_' + str(id) + '\n')
                spacer.write(spacer_dic[dr_region][id] + '\n')
                spacer_id_2_seq[dr_region + '_spacer_' + str(id)] = spacer_dic[dr_region][id]
                spacer_seq_2_dr_region[spacer_dic[dr_region][id]] = dr_region
spacer.close()

# make blast database for spacer sequences
os.system('makeblastdb -in spacers.fasta -dbtype nucl')

# targeting search
# prepare sequences for targeting search
with open('query.for.targeting.search.fasta', 'w') as q:
    for header in header2seq.keys():
        if header not in link_via_dr.keys():
            q.write('>' + header + '\n')
            q.write(header2seq[header] + '\n')
        else:
            pass
q.close()


# perform blast search
print('[{0}] [10/12] Performing BLASTn search ...'.format(strftime("%c")), end='\n', file=log, flush=True)
os.system('blastn -task blastn -db spacers.fasta -query query.for.targeting.search.fasta -out query.for.targeting.search.fasta.blastn.spacers.fasta -evalue 1e-3 -outfmt 6 '
          '-perc_identity 70 -num_threads {0}'.format(args.process))


# parse blast search results
print('[{0}] [11/12] Parsing BLASTn search results...'.format(strftime("%c")), end='\n', file=log, flush=True)
out = open('query.for.targeting.search.fasta.blastn.spacers.fasta.parsed.txt', 'w')
out.write('\t'.join(["Targeted_scaffold", "Targeted_scaffold_length", "Targeted_position", "Targeted_sequence", "Targeting_spacer", "Targeting_spacer_sequence",
                     "Targeting_spacer_origin", "Aligned_length", "Aligned_mismatches", "CRISPR_sequence", "CRISPR_sequence_Cas", "\n"]))

with open('query.for.targeting.search.fasta.blastn.spacers.fasta', 'r') as a:
    for line in a.readlines():
        line = line.strip().split('\t')
        if line[5] == '0':
            if int(line[6]) < int(line[7]):
                if spacer_id_2_seq[line[1]] in dr_region_2_read_spacer[line[1].rsplit('_spacer_', 1)[0]] and spacer_id_2_seq[line[1]] not in dr_region_2_sequence_spacer[
                    line[1].rsplit('_spacer_', 1)[0]]:
                    out.write(line[0] + '\t' + str(len(header2seq[line[0]])) + '\t' + line[6] + '-' + line[7] + '\t' + header2seq[line[0]][int(line[6]) - 1:int(line[7])] +
                                '\t' + line[1] + '\t' + spacer_id_2_seq[line[1]] + '\t' + 'read' + '\t' + line[3] + '\t' + line[4] + '\t' + line[1].rsplit('_spacer', 1)[0] +
                                '\t' + ','.join(dr_region_2_cas[line[1].rsplit('_spacer', 1)[0]][:]) + '\n')
                else:
                    out.write(line[0] + '\t' + str(len(header2seq[line[0]])) + '\t' + line[6] + '-' + line[7] + '\t' + header2seq[line[0]][int(line[6]) - 1:int(line[7])] +
                                '\t' + line[1] + '\t' + spacer_id_2_seq[line[1]] + '\t' + 'sequence,read' + '\t' + line[3] + '\t' + line[4] + '\t' + line[1].rsplit('_spacer', 1)[0] +
                                '\t' + ','.join(dr_region_2_cas[line[1].rsplit('_spacer', 1)[0]][:]) + '\n')
            else:
                if spacer_id_2_seq[line[1]] in dr_region_2_read_spacer[line[1].rsplit('_spacer_', 1)[0]] and spacer_id_2_seq[line[1]] not in dr_region_2_sequence_spacer[
                    line[1].rsplit('_spacer_', 1)[0]]:
                    out.write(line[0] + '\t' + str(len(header2seq[line[0]])) + '\t' + line[6] + '-' + line[7] + '\t' + header2seq[line[0]][int(line[7]) - 1:int(line[6])] +
                                '\t' + line[1] + '\t' + spacer_id_2_seq[line[1]] + '\t' + 'read' + '\t' + line[3] + '\t' + line[4] + '\t' + line[1].rsplit('_spacer', 1)[0] +
                                '\t' +','.join(dr_region_2_cas[line[1].rsplit('_spacer', 1)[0]][:]) + '\n')
                else:
                    out.write(line[0] + '\t' + str(len(header2seq[line[0]])) + '\t' + line[6] + '-' + line[7] + '\t' + header2seq[line[0]][int(line[7]) - 1:int(line[6])] +
                                '\t' + line[1] + '\t' + spacer_id_2_seq[line[1]] + '\t' + 'sequence,read' + '\t' + line[3] + '\t' + line[4] + '\t' + line[1].rsplit('_spacer', 1)[0] +
                                '\t' + ','.join(dr_region_2_cas[line[1].rsplit('_spacer', 1)[0]][:]) + '\n')
        else:
            pass
a.close()
out.close()


# save CRISPR-Cas system information
print('[{0}] [12/12] Saving the information of all identified CRISPR-Cas systems ...'.format(strftime("%c")), end='\n', file=log, flush=True)
with open('identified.CRISPR-Cas.system.info.txt', 'w') as info:
    title = ['Sequence', 'Repeat_region', 'Repeat_seq', 'Cas_proteins', 'Num_unique_sequence_spacers', 'Num_unique_read_spacers', '\n']
    info.write('\t'.join(title[:]))
    for scaffold in scaff2array2start_repeat.keys():
        for array in scaff2array2start_repeat[scaffold].keys():
            start = scaff2array2start_repeat[scaffold][array][0]
            end = scaff2array2start_repeat[scaffold][array][1]
            dr = scaff2array2start_repeat[scaffold][array][2]
            dr_region = scaffold + '_region_' + start + '_' + end
            if dr_region in dr_region_2_cas.keys() and dr_region in dr_region_2_sequence_spacer.keys():
                try:
                    contents = [scaffold, start + '-' + end, dr, ','.join(dr_region_2_cas[dr_region]), str(len(dr_region_2_sequence_spacer[dr_region])),
                                str(len(dr_region_2_read_spacer[dr_region])), '\n']
                    info.write('\t'.join(contents[:]))
                except KeyError:
                    print('{0} does not have spacer from read.'.format(dr_region), flush=True)
            else:
                pass
info.close()

