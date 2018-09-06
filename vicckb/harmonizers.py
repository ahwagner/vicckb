import requests
import os
import csv

BIOONTOLOGY_URL = 'http://data.bioontology.org/search'
EBI_URL = 'https://www.ebi.ac.uk/ols/api/search'
ONTOLOGIES = {
    'DOID': {'ebi': 'doid', 'bioontology': 'DOID'}
}


class DiseaseHarmonizer:

    def __init__(self, api_key=None, map_file=None, disease_ontology=None):
        if api_key is None:
            self.api_key = os.environ['BIOONTOLOGY_API_KEY']
        else:
            self.api_key = api_key
        if disease_ontology:
            assert disease_ontology in ONTOLOGIES
        self.disease_ontology = disease_ontology
        self._cache = dict()
        self._init_map(map_file)

    def _init_map(self, filename=None):
        """Subclasses may override this to preload source-specific terms.
        It may also be called with the optional filename to preload from a TSV.
        TSV is expected to be of format: source term, expanded term"""
        self._map = dict()
        if filename is None:
            return
        with open(filename, 'r', newline='') as f:
            reader = csv.reader(f, delimiter="\t")
            for line in reader:
                self._map[line[0].lower()] = line[1]

    def harmonize(self, term):
        assert isinstance(term, str)
        alias = self._map.get(term.lower(), False)
        if alias:
            term = alias
        if term.lower() in self._cache:
            return self._cache[term.lower()]
        result = self.query_ebi(term)
        if result:
            return self._prepare_result(term, result, 'ebi')
        result = self.query_bioontology(term)
        if result:
            return self._prepare_result(term, result, 'bioontology')
        self._cache[term] = None
        return None

    def _prepare_result(self, term, result, engine):
        out = {k: v for k, v in result.items() if k in ['ontology', 'term', 'id']}
        out['resultEngine'] = engine
        self._cache[term.lower()] = out
        return out

    def _submit_query(self, url, payload):
        r = requests.get(url, params=payload)
        r.raise_for_status()
        return r.json()

    def query_ebi(self, term):
        payload = {
            'q': term,
            'groupField': 'iri',
            'exact': 'on',
            'start': '0'
        }
        if self.disease_ontology and ONTOLOGIES[self.disease_ontology]['ebi']:
            payload['ontology'] = ONTOLOGIES[self.disease_ontology]['ebi']
        j = self._submit_query(EBI_URL, payload)['response']
        if j['numFound'] == 0:
            return None
        match = j['docs'][0]
        return {
            'id': match['obo_id'],
            'ontology': match['ontology_prefix'],
            'term': match['label']
        }

    def query_bioontology(self, term):
        payload = {'q': term, 'apikey': self.api_key}
        if self.disease_ontology and ONTOLOGIES[self.disease_ontology]['bioontology']:
            payload['ontologies'] = ONTOLOGIES[self.disease_ontology]['bioontology']
        j = self._submit_query(BIOONTOLOGY_URL, payload)
        if j.get('collection', False):
            match = j['collection'][0]
            prefLabel = match['prefLabel']
            ontology, oid = match['@id'].split('/')[-2:]
            if ontology == 'obo':
                ontology = oid.split('_')[0]
            return {
                'ontology': ontology,
                'term': prefLabel,
                'id': oid.replace('_', ':'),
                'matchType': match['matchType']
            }
        return None
