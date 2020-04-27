import csv
import json
import os
import random
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from numpy import median
import extract_param



# generate path for data base
DATABASE_NUMBER = str(random.randrange(1000000))
RES_DIR = 'data_base' + DATABASE_NUMBER + '/'

# random selected parameters for the data base generation
RES_DIR = 'data_base' + DATABASE_NUMBER + '/'
SUB_MODEL_DIR = 'param/' + random.choice(os.listdir('param/'))
SUB_MODELS = ['GTR', 'GTR+G', 'GTR+I+G', 'HKY', 'HKY+G', 'HKY+I+G']
SELECT_SUB_MODEL = SUB_MODELS[random.randrange(6)]
POP_SIZE = random.randrange(2, 16)
N_MORPH = random.randrange(2, 11)
D_FEATURES = random.randrange(1, N_MORPH + 1)
CONTINUOUS_SIGMA = random.randrange(151)
KMER_LEN = 5
GENETIC_SIGMA = random.randrange(61)
SIG2 = []
for n in range(N_MORPH - D_FEATURES+1):
    SIG2.append(random.uniform(0.1, 1))
BIRTH_RATE = random.uniform(0.01, 1)
DEATH_RATE = random.uniform(0.01, 0.5)
while DEATH_RATE > BIRTH_RATE:
    DEATH_RATE = random.uniform(0.009, BIRTH_RATE)
N_SPECIES = random.randrange(5, 26)
KMER = 1
MDS = 0
CF = 1
DF = 1
RUNS = 10
while KMER + MDS + CF + DF == 0:
    KMER = random.getrandbits(1)
    #MDS = random.getrandbits(1)
    CF = random.getrandbits(1)
    DF = random.getrandbits(1)

# output file
RESULT_DIR = 'results.csv'

def make_alignment(alignment_dir):
    return AlignIO.read(open(alignment_dir), 'fasta')

def get_genetic_distances(aln):
    alignment = make_alignment(aln)
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    distances = [max(max(dm)), median(dm)]
    return distances


# extract substitution parameters from files
# input: sub model path, sub model to be used, birth rate, death rate, number of species
# output: parameters used for sub model
sub_params = extract_param.extract_param(SUB_MODEL_DIR,
                                         SELECT_SUB_MODEL,
                                         BIRTH_RATE, DEATH_RATE, N_SPECIES)

# concatinate parameter string for data base generation and run generator
params = [RES_DIR, str(POP_SIZE), str(D_FEATURES), str(CONTINUOUS_SIGMA), str(KMER_LEN), str(GENETIC_SIGMA),
          str(SIG2).replace(' ', ''), SELECT_SUB_MODEL, str(sub_params).replace(' ', '')]
s = " "
test = s.join(params)
os.system('python3 data_base_gen.py ' + s.join(params))


#calculete genetic distances
distances = get_genetic_distances(RES_DIR + 'seq_gen_out/Dna.fa')

# run machine learning script
os.system(
    'python3 svm.py ' + RES_DIR + 'data_base.json 0.2 ' + str(KMER) + ' ' + str(MDS) + ' ' + str(CF) + ' ' + str(DF) + ' ' + str(RUNS))

# if result file does not exist yet, create new file with field names
if not os.path.exists(RESULT_DIR):
    with open(RESULT_DIR, 'w') as file:
        csvwriter = csv.writer(file, delimiter=',')
        csvwriter.writerow(
            ['sub model', 'base freq', 'sub/trans', 'I+G', 'n species', 'pop size', 'd traits', 'c traits', 'c trait sigma',
             'kmer size', 'distance max', 'distance median', 'genetic sigma', 'brown sigma',
             'birth rate', 'death rate', 'svm score', 'nn score', 'svm predict', 'nn predict', 'svm correct %', 'ann correct %',' p_min', 'p_max', 'kmer',
             'mds', 'cF', 'dF'])

# load results from machine learning
with open('svm.json', 'r') as svm_file:
    svm_results = json.load(svm_file)

# write used parameters and results to file
with open(RESULT_DIR, 'a') as file:
    csvwriter = csv.writer(file, delimiter=',')
    results = [SELECT_SUB_MODEL, sub_params[0], sub_params[1], sub_params[2], N_SPECIES, POP_SIZE, D_FEATURES,
               len(SIG2), CONTINUOUS_SIGMA, KMER_LEN, distances[0], distances[1], GENETIC_SIGMA, SIG2,
               BIRTH_RATE, DEATH_RATE]
    results.extend(svm_results)
    results.extend([bool(KMER), bool(MDS), bool(CF), bool(DF)])
    csvwriter.writerow(results)

# remove data base
os.system('rm -r ' + RES_DIR)