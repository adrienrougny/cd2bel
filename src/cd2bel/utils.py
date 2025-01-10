import libchebipy
import requests


def get_chebi_label_from_id(id_):
    entity = libchebipy.ChebiEntity(id_)
    name = entity.get_name()
    return name


def get_go_label_from_id(id_):
    url = f"https://api.geneontology.org/api/ontology/term/{id_}"
    response = requests.get(url)
    json = response.json()
    label = json.get("label")
    return label


def get_interpro_label_from_id(id_):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{id_}"
    response = requests.get(url)
    json = response.json()
    metadata = json.get("metadata")
    if metadata is not None:
        name_long_and_short = metadata.get("name")
        if name_long_and_short is not None:
            name_long = name_long_and_short.get("name")
            return name_long
    return None


def get_mesh_label_from_id(id_):
    url = f"https://id.nlm.nih.gov/mesh/lookup/label?resource={id_}"
    response = requests.get(url)
    json = response.json()
    if json:
        return json[0]
    return None


def get_pfam_label_from_id(id_):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{id_}"
    response = requests.get(url)
    json = response.json()
    metadata = json.get("metadata")
    if metadata is not None:
        name_long_and_short = metadata.get("name")
        if name_long_and_short is not None:
            name_long = name_long_and_short.get("name")
            return name_long
    return None


def get_reactome_label_from_id(id_):
    url = f"https://reactome.org/ContentService/data/query/{id_}"
    response = requests.get(url)
    json = response.json()
    names = json.get("name")
    if names:
        return names[0]
    return None
