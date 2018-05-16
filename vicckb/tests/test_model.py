import pytest
from collections import Counter
from vicckb.model import ViccDb, GenomicFeature, Disease
from os import environ
import networkx


CACHE_PRESENT = ViccDb.DEFAULT_CACHE.exists()
PAPER_TESTING = bool(environ.get('VICC_PAPER_TESTING'))
if PAPER_TESTING:
    SOURCES = ['molecularmatch', 'civic', 'pmkb', 'oncokb', 'jax', 'cgi']
else:
    SOURCES = ['molecularmatch', 'civic', 'pmkb', 'brca', 'oncokb', 'jax', 'cgi']


@pytest.fixture(scope="module")
def vdb():
    vdb = ViccDb(load_cache=CACHE_PRESENT)
    return vdb


@pytest.fixture(scope="module", params=SOURCES)
def sourcedb(vdb, request):
    return ViccDb(vdb.by_source(request.param))


@pytest.fixture(scope="module")
def gfa():
    return GenomicFeature(1, 1, 1, 'GRCh37', 'Feature A', 'FAKE1', alt='G')


@pytest.fixture(scope="module")
def gfb():
    return GenomicFeature(1, 1, 1, 'GRCh37', 'Feature A', 'FAKE1', alt='C')


@pytest.fixture(scope="module")
def gfa2():
    return GenomicFeature(1, 1, 1, 'GRCh37', 'Feature A', 'FAKE1', alt='G', sequence_ontology={'soid': 141})


class TestViccDb(object):

    def test_len(self, vdb):
        assert len(vdb) > 5000

    def test_select(self, vdb):
        civicdb = vdb.by_source('civic')
        assert len(civicdb) == len(vdb.select(lambda x: x['source'] == 'civic'))
        assert len(civicdb) > 1000

    def test_iter(self, vdb):
        count = 0
        for _ in vdb:
            count += 1
        assert len(vdb) == count

    def test_subtraction(self, vdb):
        civicdb = vdb.by_source('civic')
        delta = vdb - civicdb
        assert len(delta) == len(vdb) - len(civicdb)
        assert len(delta) > 5000

    def test_search_features(self, vdb):
        hits = vdb.search_by_feature(chromosome=7, start=140453136, end=140453136,
                                     reference_name='GRCh37', name='V600E', gene_symbol='BRAF')
        assert len(hits) >= 500
        associations = [hit['association'] for hit in hits]
        results = ViccDb(associations)
        assert len(results.sources) >= 5
        gf = GenomicFeature(chromosome=7, start=140453136, end=140453136,
                                     reference_name='GRCh37', name='V600E', gene_symbol='BRAF')
        hits2 = list(vdb.search_by_features([gf]))
        assert len(hits2) == len(hits)

    def test_multisearch_features(self, vdb):
        TEST_SIZE = 500
        unique_features = set()
        x = [x.features for x in vdb]
        for fset in x:
            unique_features.update(fset)
        unique_features = list(unique_features)
        hits = vdb.search_by_features(unique_features[:TEST_SIZE])
        aset = {hit['association'] for hit in hits}
        assert len(aset) > TEST_SIZE
        qset = {hit['query'] for hit in hits}
        assert len(qset) < TEST_SIZE


class TestGenomicFeatures(object):

    def test_len(self, vdb):
        for association in vdb:
            for feature in association.features:
                try:
                    assert len(feature) >= 1
                except ValueError:
                    raise ValueError("Association {} has feature {} with invalid length".format(association, feature))

    def test_equality(self, gfa, gfb):
        assert gfa == gfb

    def test_hash(self, gfa, gfb, gfa2):
        assert hash(gfa) != hash(gfb)
        assert hash(gfa) == hash(gfa2)


class TestDisease(object):

    def test_ontology(self):
        graph = Disease.DISEASE_ONTOLOGY
        assert networkx.is_directed_acyclic_graph(graph)

    def test_uniqueness(self):
        s = set()
        disease_a = Disease('DOID:1', 'DOID', 'Disease_1')
        disease_b = Disease('DOID:1', 'DOID', 'Disease 1')
        disease_c = Disease('DOID:2', 'DOID', 'Disease_2')
        s.add(disease_a)
        s.add(disease_b)
        s.add(disease_c)
        assert 2 == len(s)


class TestOncokb(object):

    def test_types(self, vdb):
        okb = vdb.by_source('oncokb')
        for x in okb:
            assert 'clinical' in x['raw']  # Only clinical results


class TestSource(object):

    def test_hash(self, sourcedb):
        c = Counter(map(lambda x: hash(x), sourcedb))   # Implicit test that all elements hash
        assert len(c) == sum(c.values())             # Tests that hashes are unique

    def test_genes(self, sourcedb):
        count = 0
        for x in sourcedb:
            if len(x.genes) == 0:
                count += 1
        assert count < 0.01 * len(sourcedb)  # less than 1% of associations lacking genes

    def test_features(self, sourcedb):
        count = 0
        for x in sourcedb:
            if len(x.features) == 0:
                count += 1
        assert count < 0.02 * len(sourcedb)  # less than 2% of associations lacking features

    def test_diseases(self, sourcedb):
        count = 0
        for x in sourcedb:
            if x.disease is None:
                count += 1
        assert count < 5

    def test_evidence_level(self, sourcedb):
        count = 0
        for x in sourcedb:
            if not x.evidence_level:
                count += 1
        assert count == 0

    # def test_drugs(self, sourcedb):
    #     count = 0
    #     for x in sourcedb:
    #         if len(x.drugs) == 0:
    #             count += 1
    #     assert count == 0