import nucpos_vars
from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from mechanical_models import probability_landscape
from energy_density import energy

def load_insertion(insertion_name, barcode) :
    """
    Loads an insertion and subsitutes the 20 barcode nucleotides with the
    given 'barcode' sequence
    """
    seqfile = '%s/%s.seq'%(nucpos_vars.sequences_dir, insertion_name)
    with open (seqfile, "r") as myfile:
        seq = myfile.read()
    return seq.replace(20*'N',barcode)

def load_Drosophila_genome() :
    """
    Load the Drosophila genome (using Biopython)
    """
    dm_genome_file = '/mnt/shared/seq/dm3R5/dmel-all-chromosome-r5.53_oneline.fasta'
    return SeqIO.index(dm_genome_file, 'fasta', alphabet=IUPAC.unambiguous_dna)

class Sequence(object) :
    """
    A Sequence objects is initialized with a string and gives the possibility to
    calculate the probability landscape, energy landscape of nucleosome
    occupancy.
    """
    
    def __init__(self,seq) :
        if not isinstance(seq,basestring) :
            raise TypeError('Init sequence with string!')
        self.seq    = seq     # sequence
        self._p     = {}      # probability density
        self._E     = {}      # energy landscape
        self.nuc    = {}      # nucleosome occupancy
        
    def p(self, order, mechanical_model, temperature) :
        if not self._p.has_key((order, mechanical_model, temperature)) :
            this_p = probability_landscape(self.seq, order, mechanical_model, temperature)
            self._p[(order,mechanical_model,temperature)] = this_p
        return self._p[(order, mechanical_model, temperature)]
    
    def E(self, order, mechanical_model, temperature) :
        if not self._E.has_key((order, mechanical_model, temperature)) :
            this_E = energy(self.p(order, mechanical_model, temperature))
            self._E[(order,mechanical_model,temperature)] = this_E
        return self._E[(order, mechanical_model, temperature)]

class InsertedSequence(Sequence) :
    """
    An InsertedSequence is a subclass of a Sequence, which takes a genome, a
    chromosome, and inserts the sequence into the chromosome at a given site.
    """
    def __init__(self,seq,genome,chromosome,cut_site,
                 left=nucpos_vars.left,right=nucpos_vars.right) :
        # set class properties
        self.chromosome = chromosome
        self.cut_site = cut_site
        self.left = left
        self.right = right
        # insert the sequence in the genome at requested position
        c = genome[chromosome]
        fullseq = c[cut_site-left:cut_site] + seq + c[cut_site:cut_site+right]
        # invoke the parent class constructor
        super(InsertedSequence, self).__init__(str(fullseq.seq))
