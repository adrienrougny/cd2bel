import momapy.celldesigner.core
import momapy.celldesigner.io.celldesigner
import momapy.io
import momapy.builder

import momapy_bel.core
import momapy_bel.io.bel

CD_NAMESPACE = "CELLDESIGNER"
PROTEIN_NAMESPACES = ["HGNC", "UNIPROT"]
ABUNDANCE_NAMESPACES = ["CHEBI"]
COMPLEX_NAMESPACES = []
LOCATION_NAMESPACES = []


def cd_to_bel(input_file_path, output_file_path):
    cd_map = momapy.io.read(input_file_path)
    cd_model = cd_map.model
    cd_annotations = cd_map.map_element_to_annotations
    bel_model = cd_model_to_bel_model(cd_model, cd_annotations)
    bel_model = momapy.builder.object_from_builder(bel_model)
    momapy.io.write(bel_model, output_file_path, "bel")
    return bel_model


def cd_model_to_bel_model(cd_model, cd_annotations):
    cd_element_to_bel_element = {}
    bel_model = momapy_bel.core.BELModelBuilder()
    for cd_species in cd_model.species:
        _ = make_and_add_bel_element_from_cd_element(
            cd_species, cd_annotations, bel_model, cd_element_to_bel_element
        )
    for cd_reaction in cd_model.reactions:
        _ = make_and_add_bel_element_from_cd_element(
            cd_reaction, cd_annotations, bel_model, cd_element_to_bel_element
        )
    bel_model = momapy.builder.object_from_builder(bel_model)
    return bel_model


def make_and_add_location_from_cd_compartment(
    cd_compartment,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element,
    super_bel_element,
):
    bel_location = cd_element_to_bel_element.get(cd_compartment)
    if bel_location is None:
        bel_location = bel_model.new_element(momapy_bel.core.Location)
        main_iri = get_main_annotation_iri_from_cd_element(
            cd_compartment, cd_annotations, namespaces=LOCATION_NAMESPACES
        )
        if main_iri is not None:
            namespace, identifier = get_namespace_and_identifier_from_cd_iri(
                main_iri
            )
        else:
            namespace = CD_NAMESPACE
            identifier = (
                cd_compartment.name if cd_compartment.name else "default"
            )
            identifier = get_normalized_identifier(identifier)
        bel_location.namespace = namespace
        bel_location.identifier = identifier
        cd_element_to_bel_element[cd_compartment] = bel_location
    super_bel_element.location = bel_location
    return bel_location


def make_and_add_abundance_from_cd_species(
    cd_species,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_abundance = bel_model.new_element(momapy_bel.core.Abundance)
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_species, cd_annotations, ABUNDANCE_NAMESPACES
    )
    if main_iri is not None:
        namespace, identifier = get_namespace_and_identifier_from_cd_iri(
            main_iri
        )
    else:
        namespace = CD_NAMESPACE
        identifier = cd_species.name
        identifier = get_normalized_identifier(identifier)
    bel_abundance.namespace = namespace
    bel_abundance.identifier = identifier
    if cd_species.compartment is not None:
        cd_compartment = cd_species.compartment
        make_and_add_location_from_cd_compartment(
            cd_compartment=cd_compartment,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_species,
            super_bel_element=bel_abundance,
        )
    bel_abundance = momapy.builder.object_from_builder(bel_abundance)
    bel_model.statements.add(bel_abundance)
    cd_element_to_bel_element[cd_species] = bel_abundance
    return bel_abundance


def make_and_add_modification_from_cd_modification(
    cd_modification,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_modification = bel_model.new_element(
        momapy_bel.core.ProteinModification
    )
    bel_modification.identifier = cd_modification.state.value
    if (
        cd_modification.residue is not None
        and cd_modification.residue.name is not None
    ):
        cd_residue = cd_modification.residue.name
        if not cd_residue[0].isnumeric():
            bel_modification.amino_acid = cd_residue[0]
            bel_modification.residue = cd_residue[1:]
        else:
            bel_modification.residue = cd_residue
    bel_modification = momapy.builder.object_from_builder(bel_modification)
    super_bel_element.modifications.append(bel_modification)
    return bel_modification


def make_and_add_protein_abundance_from_cd_species(
    cd_species,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_abundance = bel_model.new_element(momapy_bel.core.ProteinAbundance)
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_species, cd_annotations, PROTEIN_NAMESPACES
    )
    if main_iri is not None:
        namespace, identifier = get_namespace_and_identifier_from_cd_iri(
            main_iri
        )
    else:
        namespace = CD_NAMESPACE
        identifier = cd_species.name
        identifier = get_normalized_identifier(identifier)
    bel_abundance.namespace = namespace
    bel_abundance.identifier = identifier
    if cd_species.compartment is not None:
        cd_compartment = cd_species.compartment
        make_and_add_location_from_cd_compartment(
            cd_compartment=cd_compartment,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_species,
            super_bel_element=bel_abundance,
        )
    for cd_modification in cd_species.modifications:
        if (
            cd_modification.state is not None
            and cd_modification.state.value != "empty"
        ):  # to delete second part
            make_and_add_modification_from_cd_modification(
                cd_modification=cd_modification,
                cd_annotations=cd_annotations,
                bel_model=bel_model,
                cd_element_to_bel_element=cd_element_to_bel_element,
                super_cd_element=cd_species,
                super_bel_element=bel_abundance,
            )
    bel_abundance = momapy.builder.object_from_builder(bel_abundance)
    bel_model.statements.add(bel_abundance)
    cd_element_to_bel_element[cd_species] = bel_abundance
    return bel_abundance


def make_and_add_complex_abundance_from_cd_species(
    cd_species,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_abundance = bel_model.new_element(momapy_bel.core.ComplexAbundance)
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_species, cd_annotations, COMPLEX_NAMESPACES
    )
    if main_iri is not None:
        namespace, identifier = get_namespace_and_identifier_from_cd_iri(
            main_iri
        )
    else:
        namespace = CD_NAMESPACE
        identifier = cd_species.name
        identifier = get_normalized_identifier(identifier)
    bel_abundance.namespace = namespace
    bel_abundance.identifier = identifier
    if cd_species.compartment is not None:
        cd_compartment = cd_species.compartment
        make_and_add_location_from_cd_compartment(
            cd_compartment=cd_compartment,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_species,
            super_bel_element=bel_abundance,
        )
    for cd_subunit in cd_species.subunits:
        make_and_add_bel_element_from_cd_element(
            cd_element=cd_subunit,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_species,
            super_bel_element=bel_abundance,
        )
    bel_abundance = momapy.builder.object_from_builder(bel_abundance)
    bel_model.statements.add(bel_abundance)
    cd_element_to_bel_element[cd_species] = bel_abundance
    return bel_abundance


def make_and_add_included_abundance_from_cd_subunit(
    cd_subunit,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element,
    super_bel_element,
):
    bel_abundance = bel_model.new_element(momapy_bel.core.Abundance)
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_subunit, cd_annotations, ABUNDANCE_NAMESPACES
    )
    if main_iri is not None:
        namespace, identifier = get_namespace_and_identifier_from_cd_iri(
            main_iri
        )
    else:
        namespace = CD_NAMESPACE
        identifier = cd_subunit.name
        identifier = get_normalized_identifier(identifier)
    bel_abundance.namespace = namespace
    bel_abundance.identifier = identifier
    if cd_subunit.compartment is not None:
        cd_compartment = cd_subunit.compartment
        make_and_add_location_from_cd_compartment(
            cd_compartment=cd_compartment,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_subunit,
            super_bel_element=bel_abundance,
        )
    bel_abundance = momapy.builder.object_from_builder(bel_abundance)
    super_bel_element.members.add(bel_abundance)
    return bel_abundance


def make_and_add_reaction_from_cd_reaction(
    cd_reaction,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    reactants = []
    products = []
    for cd_reactant in cd_reaction.reactants:
        cd_reactant_species = cd_reactant.referred_species
        bel_reactant = cd_element_to_bel_element.get(cd_reactant_species)
        if bel_reactant is None:
            bel_reactant = make_and_add_bel_element_from_cd_element(
                cd_reactant_species,
                cd_annotations,
                bel_model,
                cd_element_to_bel_element,
            )
        if bel_reactant is not None:  # to delete
            reactants.append(bel_reactant)
    for cd_product in cd_reaction.products:
        cd_product_species = cd_product.referred_species
        bel_product = cd_element_to_bel_element.get(cd_product_species)
        if bel_product is None:
            bel_product = make_and_add_bel_element_from_cd_element(
                cd_product_species,
                cd_annotations,
                bel_model,
                cd_element_to_bel_element,
            )
        if bel_product is not None:  # to delete
            products.append(bel_product)
    if reactants and products:  # to delete
        bel_reaction = momapy_bel.core.Reaction(
            reactants=frozenset(reactants), products=frozenset(products)
        )
        bel_model.statements.add(bel_reaction)
        cd_element_to_bel_element[cd_reaction] = bel_reaction
    else:
        bel_reaction = None
    return bel_reaction


def make_and_add_modulatory_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_model=bel_model,
        bel_model_function_name="regulates",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect=None,
    )


def make_and_add_catalytic_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_model=bel_model,
        bel_model_function_name="add_directly_increases",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect={"namespace": "GO", "identifier": "0003824"},
    )


def make_and_add_inhibitory_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_model=bel_model,
        bel_model_function_name="add_decreases",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect=None,
    )


def make_and_add_stimulatory_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element,
):
    generic_make_and_add_activity_from_cd_modifier(
        cd_modifier=cd_modifier,
        cd_annotations=cd_annotations,
        bel_model=bel_model,
        bel_model_function_name="add_directly_increases",
        cd_element_to_bel_element=cd_element_to_bel_element,
        super_cd_element=super_cd_element,
        effect=None,
    )


def generic_make_and_add_activity_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_model,
    bel_model_function_name,
    cd_element_to_bel_element,
    super_cd_element,
    effect=None,
):
    bel_reaction = cd_element_to_bel_element[super_cd_element]
    cd_species = cd_modifier.referred_species
    bel_abundance = cd_element_to_bel_element.get(cd_species)
    if bel_abundance is None:
        bel_abundance = make_and_add_bel_element_from_cd_element(
            cd_species, cd_annotations, bel_model, cd_element_to_bel_element
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
        for cd_modification in cd_species.modifications:
            print(cd_modification)
        if effect is not None:
            source_modifier["effect"] = effect
        bel_model_function = getattr(bel_model, bel_model_function_name)
        bel_model_function(
            source=bel_abundance,
            target=bel_reaction,
            evidence="test evidence",
            citation="12345678",
            source_modifier=source_modifier,
        )


def make_and_add_bel_element_from_cd_element(
    cd_element,
    cd_annotations,
    bel_model,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    make_and_add_function = get_make_and_add_function_from_cd_element(
        cd_element, super_cd_element
    )
    if make_and_add_function is not None:
        bel_term = make_and_add_function(
            cd_element,
            cd_annotations,
            bel_model,
            cd_element_to_bel_element,
            super_cd_element,
            super_bel_element,
        )
    else:
        bel_term = None
        # print(f"no function for element of type {type(cd_element)}")
    return bel_term


def get_make_and_add_function_from_cd_element(cd_element, super_cd_element):
    if super_cd_element is not None:
        key = (type(cd_element), type(super_cd_element))
    else:
        key = type(cd_element)
    make_and_add_function = CD_ELEMENT_TYPE_TO_MAKE_AND_ADD_FUNC.get(key)
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
    identifier = get_normalized_identifier(identifier)
    return namespace, identifier


def get_normalized_identifier(cd_identifier):
    identifier = cd_identifier.replace("GO%3A", "")
    identifier = identifier.replace("CHEBI%3A", "")
    if not identifier.isalnum():
        identifier = f'"{identifier}"'
    return identifier


def get_normalized_namespace(namespace):
    if namespace in CD_NAMESPACE_TO_NORMALIZED_NAMESPACE:
        return CD_NAMESPACE_TO_NORMALIZED_NAMESPACE[namespace]
    else:
        return namespace


CD_ELEMENT_TYPE_TO_MAKE_AND_ADD_FUNC = {
    momapy.celldesigner.core.SimpleMolecule: make_and_add_abundance_from_cd_species,
    momapy.celldesigner.core.Complex: make_and_add_complex_abundance_from_cd_species,
    momapy.celldesigner.core.Ion: make_and_add_abundance_from_cd_species,
    momapy.celldesigner.core.GenericProtein: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.Receptor: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.IonChannel: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.TruncatedProtein: make_and_add_protein_abundance_from_cd_species,
    momapy.celldesigner.core.StateTransition: make_and_add_reaction_from_cd_reaction,
    momapy.celldesigner.core.Catalyzer: make_and_add_catalytic_activity_from_cd_modifier,
    momapy.celldesigner.core.Inhibitor: make_and_add_inhibitory_activity_from_cd_modifier,
    momapy.celldesigner.core.PhysicalStimulator: make_and_add_stimulatory_activity_from_cd_modifier,
    momapy.celldesigner.core.Modulator: make_and_add_modulatory_activity_from_cd_modifier,
    (
        momapy.celldesigner.core.SimpleMolecule,
        momapy.celldesigner.core.Complex,
    ): make_and_add_included_abundance_from_cd_subunit,
}

CD_NAMESPACE_TO_NORMALIZED_NAMESPACE = {
    "hgnc.symbol": "HGNC",
    "obo.go": "GO",
    "obo.chebi": "CHEBI",
    "uniprot": "UNIPROT",
}
