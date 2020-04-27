import ast
import json
import os
import csv
import random
import sys
from Bio import AlignIO
from Bio import SeqIO
from copy import deepcopy
from kmer_counter import kmer_counter
from MDS import mds

# database parameters

RES_DIR = sys.argv[1]
POP_SIZE = int(sys.argv[2])
D_FEATURES = int(sys.argv[3])
CONTINUOUS_SIGMA = float(sys.argv[4])
KMER_LEN = int(sys.argv[5])
GENETIC_SIGMA = int(sys.argv[6])
SIG2 = ast.literal_eval(sys.argv[7])
C_FEATURES = len(SIG2)
SUB_MOD = sys.argv[8]
SUB_PARAMS = sys.argv[9]
SUB_PARAMS = ast.literal_eval(SUB_PARAMS)

# directories for database structure
SEQ_GEN_OUT = RES_DIR + 'seq_gen_out/'
GENOME_DIR = RES_DIR + 'species_genomes/'
DB_DIR = RES_DIR + 'alf/results/alf/DB/'
TREE_DIR = RES_DIR + 'RealTree.nwk'
INDIVIDUALS_GENOMES_DIR = RES_DIR + 'individuals_genomes/'
CONTINUOUS_DIR = RES_DIR + 'continious_features/'
DISCRETE_DIR = RES_DIR + 'discrete_features/'
ALIGN_DIR = RES_DIR + 'alignment/'


# adjustment for easier float calculation (calculating only with whole numbers by multiplying everything with a multiple of 10)
# input: to be modified number
# output: multiple of 10 by which to be modified
def calc_modify(x):
    s = str(x)
    if int(float(x)) > 0 and float(x) / int(float(x)) == 1:
        return 1

    i = 0
    value = 1
    while i < len(s) - s.index('.') - 1:
        value *= 10
        i += 1
    return value


# return a random base for individual generation with a propability matrix
# input: base to be exchanged
# output: new base
def get_base(base):
    ref_base = base
    csv_file = open('base_ex_rate.csv', 'r')
    matrix = csv.DictReader(csv_file, delimiter=' ')
    random_number = random.randrange(100)
    return_base = ['A', 'T', 'G', 'C']
    i = 0
    prob = 0
    for row in matrix:
        prob += float(row[ref_base]) * 100
        if random_number <= prob:
            return return_base[i]
        i += 1
    csv_file.close()


# check if new generated individual sequnce is not acceptable i. e. not to diverse in regard to the original sequence
# input: original sequence, new generated sequence
# output: boolean whether accepted or not
def is_accepted(sequence, new_sequnece):
    i = 0
    diff = 0
    while i < len(sequence):
        if sequence[i] != new_sequnece[i]:
            diff += 1
        i += 1

    if diff / len(sequence) * 100 <= GENETIC_SIGMA:
        return True
    return False


# create new directory for database
# input: directory path to be created
def make_dir(dir):
    if os.path.isdir(dir) is not True:
        os.system('mkdir ' + dir)


# concatenate sequence loci to one continious sequence
# input: input directory of sequence loci, output directory of concatenated sequence
def make_genome(input, output):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input, 'fasta'))

    output = open(output, 'w')
    id = input.split('/')
    output.write('>' + id[-1] + '\n')

    for key in fasta_dict:
        output.write(str(fasta_dict[key].seq))

    output.close()
    pass

# generate sequences for individuals
# input: original sequence, output path, population size
# output: new sequence file
def gen_dna_individuals(input, output, pop_size):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input, 'fasta'))
    output = open(output, 'w')

    for key in fasta_dict:
        i = 0
        while i < pop_size:
            sequence = list(fasta_dict[key].seq)
            n_sub = len(sequence)
            while n_sub >= len(sequence):
                n_sub = abs(int(random.normalvariate(0, len(sequence) * GENETIC_SIGMA / 100)))
            ex_positions = random.sample(range(len(sequence)), n_sub)
            for n in range(len(ex_positions)):
                sequence[ex_positions[n]] = get_base(sequence[ex_positions[n]])


            output.write('>' + key + '_individual' + str(i + 1) + '\n' + ''.join(sequence) + '\n')
            i += 1

    output.close()
    pass


# read sequnce from fasta file
# input: directory of fasta file, id of individual
# output: requested sequence
def get_Genome(genome_dir, individual):
    fasta_data = list(SeqIO.parse(genome_dir, 'fasta'))
    sequence = str(fasta_data[individual].seq)
    return sequence


# read continious features from json file
# input: json file directory, species, individual
# output: list of continuous features
def get_c_features(c_dir, species, individual):
    features = []
    for file in os.listdir(c_dir):
        if file.endswith('.json'):
            data = json.load(open(c_dir + file))
            features.append(data[species][individual])
    return features


# pars csv file to dictionary
# input: csv directory
# output: dictionary
def pars_csv(input):
    data = {}
    file = open(input, 'r')
    next(file)
    for row in file:
        key = row.split(',')[0].replace('"', '')
        features = []
        data.update({key: features})
        data[key].append(float(row.split(',')[1].replace('\n', '').replace('"', '')))
    file.close()
    return data


# read discrete features from csv file
# input: discrete feature csv directory, species
# output: list of discrete features
def get_d_features(d_dir, species):
    features = []
    for dir in os.listdir(d_dir):
        for file in os.listdir(d_dir + '/' + dir):
            if file.endswith('.csv'):
                data = pars_csv(d_dir + '/' + dir + '/' + file)
                features.extend(data[species])
    return features


# save database as json file
# input: list of species names, kmer data, mds features
def make_json(species_names, kmer_data, mds_features):
    json_file = open(RES_DIR + 'data_base.json', 'w')
    species_count = 0
    data = {'kmer_list': kmer_data['kmer_list']}
    for species in species_names:
        individuals = []
        i = 0
        while i < POP_SIZE:
            genome_dir = INDIVIDUALS_GENOMES_DIR + 'pop_' + species + '_dna.fa'
            features = {'genome': get_Genome(genome_dir, i),
                        'c_features': get_c_features(CONTINUOUS_DIR, species, i),
                        'd_features': get_d_features(DISCRETE_DIR, species),
                        'kmer_profile': kmer_data[species + '_dna.fa_individual' + str(i + 1)],
                        'mds_features': list(mds_features[species_count * POP_SIZE + i])}
            individuals.append(features)
            i += 1

        data.update({species: individuals})
        species_count += 1

    json.dump(data, json_file)
    json_file.close()
    pass


# randomise continuous features to create individuals
# input: continuous feature csv directory, individual json output directory, population size, range of randomisation
def gen_con_individuals(input, output, pop_size, sta):
    data = pars_csv(input)
    for key in data:
        i = 0
        while i < pop_size - 1:
            sta_sig = data[key][0] * sta / 100
            data[key].append(data[key][0] + random.normalvariate(0, sta_sig))
            i += 1
    output_json = open(output, 'w')
    json.dump(data, output_json)
    pass


# read alignment fasta and return alignment opbject
# input: alignment fasta directory
# output: alignment object
def make_alignment(alignment_dir):
    return AlignIO.read(open(alignment_dir), 'fasta')


# parses list to string
# input: list
# output: string
def parse_to_string(string_list):
    if len(string_list) == 0:
        return ''
    string = str(string_list[0])
    for i in range(1, len(string_list)):
        string += ',' + str(string_list[i])
    return string


# splits fasta file to muliple files, one file per record and returns list of IDs
# input: path to fasta file, output path
# output: id list
def split_fasta(fasta_path, out_path):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))
    species_names = []
    for key in fasta_dict:
        species_names.append(fasta_dict[key].id)
        output = open(out_path + fasta_dict[key].id + '_dna.fa', 'w')
        sequence = str(fasta_dict[key].seq)
        output.write('>' + fasta_dict[key].id + '_dna.fa' + '\n' + sequence.upper())
        output.close()
    return species_names




# adjust alf parameter file to substitution model and run alf simulation
make_dir(RES_DIR)
test = parse_to_string(SUB_PARAMS[0])

species_names = []


# either use alf or seq-gen for sequence generarion
make_dir(SEQ_GEN_OUT)
make_dir(GENOME_DIR)

os.system('Rscript tree_sequence_simulator.R ' + SUB_MOD + ' ' + parse_to_string(SUB_PARAMS[0]) + ' ' +
          parse_to_string(SUB_PARAMS[1]) + ' ' + parse_to_string(SUB_PARAMS[2]) + ' ' + str(SUB_PARAMS[3]) + ' ' +
          str(SUB_PARAMS[4]) + ' ' + str(SUB_PARAMS[5]) + ' ' + SEQ_GEN_OUT + ' ' + TREE_DIR)

species_names = split_fasta(SEQ_GEN_OUT + 'Dna.fa', GENOME_DIR)

print('sequences done')



# generate individuals for every species
make_dir(INDIVIDUALS_GENOMES_DIR)

for file in os.listdir(GENOME_DIR):
    gen_dna_individuals(GENOME_DIR + file, INDIVIDUALS_GENOMES_DIR + 'pop_' + file, POP_SIZE)

print('individuals done')

# make alignment of individuals and generate MDS features
make_dir(ALIGN_DIR)

all_genomes = open(ALIGN_DIR + 'all_genomes.fa', 'w')
for file in os.listdir(INDIVIDUALS_GENOMES_DIR):
    handle = open(INDIVIDUALS_GENOMES_DIR + file)
    for line in handle:
        all_genomes.write(line)
    handle.close()
all_genomes.close()

alignment = make_alignment(ALIGN_DIR + 'all_genomes.fa')
mds_features = mds(alignment)

print('mds done')

# generate continuous features with R script
make_dir(CONTINUOUS_DIR)

i = 0
while i < C_FEATURES:
    os.system('Rscript fastBM.R ' + TREE_DIR + ' ' + CONTINUOUS_DIR + str(i + 1) + 'c_features.csv ' + str(SIG2[i]))
    i += 1

for file in os.listdir(CONTINUOUS_DIR):
    if file.endswith('.csv'):
        gen_con_individuals(CONTINUOUS_DIR + file, CONTINUOUS_DIR + 'pop_' + file.replace('.csv', '.json'), + POP_SIZE,
                            CONTINUOUS_SIGMA)
print('c_features done')

# generate discrete features with R script
make_dir(DISCRETE_DIR)

i = 0
while i < D_FEATURES:
    make_dir(DISCRETE_DIR + 'feature' + str(i + 1) + '/')
    os.system(
        'Rscript sim_history.R ' + TREE_DIR + ' discrete_matrix.csv ' + DISCRETE_DIR + 'feature' + str(i + 1) + '/')
    i += 1

print('d_features done')

# genearte kmer features
kmer_data = kmer_counter(INDIVIDUALS_GENOMES_DIR, KMER_LEN)
print('kmer done')

# save database as json file
make_json(species_names, kmer_data, mds_features)
