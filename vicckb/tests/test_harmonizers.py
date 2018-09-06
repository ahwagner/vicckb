import pytest
from vicckb.definitions import DATA_ROOT
from vicckb.harmonizers import DiseaseHarmonizer


@pytest.fixture(scope="module")
def adh():
    """Aliased Disease Harmonizer"""
    return DiseaseHarmonizer(map_file=(DATA_ROOT / 'disease_alias.tsv'),
                             disease_ontology='DOID')


class TestDiseaseHarmonizer(object):

    def test_adh_init(self, adh):
        assert adh

    def test_adh_aliases(self, adh):
        assert 325 == len(adh._map)

    def test_search(self, adh):
        result = adh.harmonize('breast cancer')
        expected = {
            'ontology': 'DOID',
            'term': 'breast cancer',
            'id': 'DOID:1612',
            'resultEngine': 'ebi'
        }
        assert expected == result
