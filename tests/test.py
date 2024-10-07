import cd2bel.core
import cd2bel.neo4j

import momapy.celldesigner.core
import momapy.celldesigner.io.celldesigner
import momapy.io

import pybel
import py2neo

INPUT_FILE_PATH = "/home/rougny/research/commute/commute_dm_develop/build/data/covid_dm/celldesigner/JNK_pathway.xml"
OUTPUT_FILE_PATH = "test.bel"

# cd_map = momapy.io.read(INPUT_FILE_PATH)
# cd_model = cd_map.model
# cd_annotations = cd_map.map_element_to_annotations
# bel_graph = cd2bel.core.cd_model_to_bel_graph(
#     cd_model, cd_annotations, name="test"
# )
# pybel.dump(bel_graph, OUTPUT_FILE_PATH)
bel_graph = pybel.load(OUTPUT_FILE_PATH, no_identifier_validation=True)
print(bel_graph.nodes)
print()
for edge in bel_graph.edges:
    print(edge)
