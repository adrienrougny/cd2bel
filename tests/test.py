import cd2bel.core
import cd2bel.neo4j

import momapy.celldesigner.core
import momapy.celldesigner.io.celldesigner
import momapy.io
import momapy_bel.io.bel


INPUT_FILE_PATH = "/home/rougny/research/commute/commute_dm_develop/build/data/covid_dm/celldesigner/HMOX1_pathway.xml"
OUTPUT_FILE_PATH = "test.bel"

cd_map = momapy.io.read(INPUT_FILE_PATH)
cd_model = cd_map.model
cd_annotations = cd_map.map_element_to_annotations
bel_model, bel_annotations = cd2bel.core.cd_model_to_bel_model(
    cd_model, cd_annotations
)
momapy.io.write(
    bel_model, OUTPUT_FILE_PATH, "bel", annotations=bel_annotations
)
