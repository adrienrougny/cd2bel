import momapy.celldesigner.core
import momapy.celldesigner.io.celldesigner
import momapy.io
import pybel

PROTEIN_NAMESPACES = ["HGNC"]


def cd_to_bel(input_file_path, output_file_path, name=None):
    cd_map = momapy.io.read(input_file_path)
    cd_model = cd_map.model
    cd_annotations = cd_map.map_element_to_annotations
    bel_graph = cd_model_to_bel_graph(cd_model, cd_annotations, name=name)
    pybel.dump(bel_graph, output_file_path)
    return bel_graph


def cd_model_to_bel_graph(cd_model, cd_annotations, name):
    cd_element_to_bel_element = {}
    bel_graph = pybel.BELGraph(
        name=name, description="Graph generated with the cd2bel tool"
    )
    for cd_species in cd_model.species:
        _ = make_and_add_bel_element_from_cd_element(
            cd_species, cd_annotations, bel_graph, cd_element_to_bel_element
        )
    for cd_reaction in cd_model.reactions:
        _ = make_and_add_bel_element_from_cd_element(
            cd_reaction, cd_annotations, bel_graph, cd_element_to_bel_element
        )
    bel_graph.graph[pybel.constants.GRAPH_NAMESPACE_URL] = {
        "GO": "https://raw.githubusercontent.com/pharmacome/conso/d9d270e11aac480542c412d4222983a5f042b8ae/external/go-names.belns"
    }
    return bel_graph


def make_and_add_protein_abundance_from_cd_species(
    cd_species,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element=None,
):
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_species, cd_annotations, PROTEIN_NAMESPACES
    )
    if main_iri is not None:
        namespace, identifier = get_namespace_and_identifier_from_cd_iri(
            main_iri
        )
        bel_abundance = pybel.dsl.Protein(
            namespace=namespace, identifier=identifier
        )
        bel_graph.add_node_from_data(bel_abundance)
        return bel_abundance
    return None


def make_and_add_reaction_from_cd_reaction(
    cd_reaction,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element=None,
):
    reactants = []
    products = []
    for cd_reactant in cd_reaction.reactants:
        bel_reactant = make_and_add_bel_element_from_cd_element(
            cd_reactant.referred_species,
            cd_annotations,
            bel_graph,
            cd_element_to_bel_element,
        )
        if bel_reactant is not None:  # to delete
            reactants.append(bel_reactant)
    for cd_product in cd_reaction.products:
        bel_product = make_and_add_bel_element_from_cd_element(
            cd_product.referred_species,
            cd_annotations,
            bel_graph,
            cd_element_to_bel_element,
        )
        if bel_product is not None:  # to delete
            products.append(bel_product)
    if reactants and products:
        bel_reaction = pybel.dsl.Reaction(
            reactants=reactants, products=products
        )
        bel_graph.add_node_from_data(bel_reaction)
        cd_element_to_bel_element[cd_reaction] = bel_reaction
        for cd_modifier in cd_reaction.modifiers:
            _ = make_and_add_bel_element_from_cd_element(
                cd_modifier,
                cd_annotations,
                bel_graph,
                cd_element_to_bel_element,
                cd_reaction,
            )
    else:
        bel_reaction = None
    return bel_reaction


def make_and_add_modulatory_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_graph=bel_graph,
        bel_graph_function_name="regulates",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect=None,
    )


def make_and_add_catalytic_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_graph=bel_graph,
        bel_graph_function_name="add_directly_increases",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect={"namespace": "GO", "identifier": "0003824"},
    )


def make_and_add_inhibitory_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_graph=bel_graph,
        bel_graph_function_name="add_decreases",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect=None,
    )


def make_and_add_stimulatory_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_graph=bel_graph,
        bel_graph_function_name="add_directly_increases",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect=None,
    )


def generic_make_and_add_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_graph,
    bel_graph_function_name,
    cd_element_to_bel_element,
    super_cd_element,
    effect=None,
):
    bel_reaction = cd_element_to_bel_element[super_cd_element]
    cd_species = cd_modifier.referred_species
    bel_abundance = cd_element_to_bel_element.get(cd_species)
    if bel_abundance is None:
        bel_abundance = make_and_add_bel_element_from_cd_element(
            cd_species, cd_annotations, bel_graph, cd_element_to_bel_element
        )
    if bel_abundance is not None:  # to delete
        source_modifier = {
            "modifier": "Activity",
        }
        if cd_species.compartment is not None:
            cd_compartment = cd_species.compartment
            cd_compartment_iri = get_main_annotation_iri_from_cd_element(
                cd_compartment, cd_annotations, namespaces=[]
            )
            if cd_compartment_iri is not None:
                compartment_namespace, compartment_identifier = (
                    get_namespace_and_identifier_from_cd_iri(
                        cd_compartment_iri
                    )
                )
                source_modifier["location"] = {
                    "namespace": compartment_namespace,
                    "name": compartment_identifier,
                }
            else:
                pass
        print(cd_species.species_template)
        for cd_modification in cd_species.modifications:
            print(cd_modification)
        if effect is not None:
            source_modifier["effect"] = effect
        bel_graph_function = getattr(bel_graph, bel_graph_function_name)
        bel_graph_function(
            source=bel_abundance,
            target=bel_reaction,
            evidence="test evidence",
            citation="12345678",
            source_modifier=source_modifier,
        )


def make_and_add_bel_element_from_cd_element(
    cd_element,
    cd_annotations,
    bel_graph,
    cd_element_to_bel_element,
    super_cd_element=None,
):
    make_and_add_function = get_make_and_add_function_from_cd_element(
        cd_element
    )
    if make_and_add_function is not None:
        bel_term = make_and_add_function(
            cd_element,
            cd_annotations,
            bel_graph,
            cd_element_to_bel_element,
            super_cd_element,
        )
    else:
        bel_term = None
        # print(f"no function for element of type {type(cd_element)}")
    return bel_term


def get_make_and_add_function_from_cd_element(cd_element):
    make_and_add_function = CD_ELEMENT_TYPE_TO_MAKE_AND_ADD_FUNC.get(
        type(cd_element)
    )
    return make_and_add_function


def get_main_annotation_iri_from_cd_element(
    cd_element, cd_annotations, namespaces
):
    if cd_element in cd_annotations:
        for namespace in namespaces:
            for cd_annotation in cd_annotations[cd_element]:
                cd_resources = cd_annotation.resources
                for cd_resource in cd_resources:
                    cd_namespace, _ = get_namespace_and_identifier_from_cd_iri(
                        cd_resource
                    )
                    if namespace == cd_namespace:
                        return cd_resource
        for cd_annotation in cd_annotations[
            cd_element
        ]:  # we return the first one
            for cd_resource in cd_annotation.resources:
                return cd_resource
    return None


def get_namespace_and_identifier_from_cd_iri(cd_iri):
    iri_split = cd_iri.split(":")
    namespace = iri_split[2]
    namespace = get_normalized_namespace(namespace)
    identifier = iri_split[3]
    identifier = identifier.replace("GO%3A", "")
    return namespace, identifier


def get_normalized_namespace(namespace):
    if namespace in CD_NAMESPACE_TO_NORMALIZED_NAMESPACE:
        return CD_NAMESPACE_TO_NORMALIZED_NAMESPACE[namespace]
    else:
        return namespace


CD_ELEMENT_TYPE_TO_MAKE_AND_ADD_FUNC = {
    momapy.celldesigner.core.GenericProtein: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.Receptor: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.IonChannel: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.TruncatedProtein: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.StateTransition: make_and_add_reaction_from_cd_reaction,
    momapy.celldesigner.core.Catalyzer: make_and_add_catalytic_activity_from_cd_modifier,
    momapy.celldesigner.core.Inhibitor: make_and_add_inhibitory_activity_from_cd_modifier,
    momapy.celldesigner.core.PhysicalStimulator: make_and_add_stimulatory_activity_from_cd_modifier,
    momapy.celldesigner.core.Modulator: make_and_add_modulatory_activity_from_cd_modifier,
}

CD_NAMESPACE_TO_NORMALIZED_NAMESPACE = {
    "hgnc.symbol": "HGNC",
    "obo.go": "GO",
}
