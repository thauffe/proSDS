import os
from sklearn import manifold

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# calulate distanc matrix of alignment and perform mds
# input: path to alignment
# output: mds result
def mds(aln):
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    print('distance matric done')
    man = manifold.MDS(n_components=5, metric=True, n_init=4, max_iter=300, verbose=0, eps=0.001, n_jobs=1,
                       random_state=None, dissimilarity='precomputed')
    return man.fit_transform(dm)
