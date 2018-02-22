import os

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import utils


class TestUtils(IntegronTest):

    def test_read_single_dna_fasta(self):
        replicon_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        replicon = utils.read_single_dna_fasta(replicon_path)
        self.assertEqual(replicon.id, 'ACBA.007.P01_13')
        self.assertEqual(replicon.name, replicon.name)
        self.assertEqual(len(replicon), 20301)
        self.assertTrue(replicon.seq.startswith('TGCTGCTTGGATGCCCGAGGCATAGACTGTACAAAAAAACAGTCATAACAAGCCATGAAA'))
        self.assertTrue(replicon.seq.endswith('CGACCCACGGCGTAACGCGCT'))


