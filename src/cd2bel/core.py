import collections

import momapy_bel.core
import momapy_bel.io.bel

import momapy.celldesigner.core
import momapy.celldesigner.io.celldesigner
import momapy.io
import momapy.builder

import cd2bel.utils

PROTEIN_NAMESPACES = [
    "hgnc.symbol",
    "uniprot",
    "hgnc",
    "ncbigene",
    "ensembl",
    "refseq",
    "obo.chebi",
    "ec-code",
    "ncbiprotein",
    "interpro",
    "mesh",
    "obo.go",
    "brenda",
]
SIMPLE_CHEMICAL_NAMESPACES = ["obo.chebi"]
DRUG_NAMESPACES = ["obo.chebi"]
GENE_NAMESPACES = [
    "hgnc.symbol",
    "uniprot",
    "hgnc",
    "ncbigene",
    "ensembl",
    "refseq",
    "obo.chebi",
    "ec-code",
    "ncbiprotein",
    "interpro",
    "mesh",
    "obo.go",
    "brenda",
]
RNA_NAMESPACES = [
    "hgnc.symbol",
    "uniprot",
    "hgnc",
    "ncbigene",
    "ensembl",
    "refseq",
    "obo.chebi",
    "ec-code",
    "ncbiprotein",
    "interpro",
    "mesh",
    "obo.go",
    "brenda",
]
MICRORNA_NAMESPACES = [
    "hgnc.symbol",
    "uniprot",
    "hgnc",
    "ncbigene",
    "ensembl",
    "refseq",
    "obo.chebi",
    "ec-code",
    "ncbiprotein",
    "interpro",
    "mesh",
    "obo.go",
    "brenda",
]
COMPLEX_NAMESPACES = ["obo.go"]
PHENOTYPE_NAMESPACES = ["obo.go", "mesh", "omim", "wikipathways"]
COMPARTMENT_NAMESPACES = ["obo.go", "mesh"]
PUBLICATION_NAMESPACES = ["pubmed", "doi"]
CD_NAMESPACE = "CELLDESIGNER"

CD_CLASS_TO_BEL_CLASS = {
    momapy.celldesigner.core.Degraded: momapy_bel.core.Abundance,
    momapy.celldesigner.core.SimpleMolecule: momapy_bel.core.Abundance,
    momapy.celldesigner.core.Unknown: momapy_bel.core.Abundance,
    momapy.celldesigner.core.Drug: momapy_bel.core.Abundance,
    momapy.celldesigner.core.Ion: momapy_bel.core.Abundance,
    momapy.celldesigner.core.Gene: momapy_bel.core.GeneAbundance,
    momapy.celldesigner.core.RNA: momapy_bel.core.RNAAbundance,
    momapy.celldesigner.core.AntisenseRNA: momapy_bel.core.RNAAbundance,
    momapy.celldesigner.core.Complex: momapy_bel.core.ComplexAbundance,
    momapy.celldesigner.core.GenericProtein: momapy_bel.core.ProteinAbundance,
    momapy.celldesigner.core.Receptor: momapy_bel.core.ProteinAbundance,
    momapy.celldesigner.core.IonChannel: momapy_bel.core.ProteinAbundance,
    momapy.celldesigner.core.TruncatedProtein: momapy_bel.core.ProteinAbundance,
    momapy.celldesigner.core.Phenotype: momapy_bel.core.BiologicalProcess,
    momapy.celldesigner.core.StateTransition: momapy_bel.core.Reaction,
    momapy.celldesigner.core.KnownTransitionOmitted: momapy_bel.core.Reaction,
    momapy.celldesigner.core.UnknownTransition: momapy_bel.core.Reaction,
    momapy.celldesigner.core.HeterodimerAssociation: momapy_bel.core.Reaction,
    momapy.celldesigner.core.Dissociation: momapy_bel.core.Reaction,
    momapy.celldesigner.core.Truncation: momapy_bel.core.Reaction,
    momapy.celldesigner.core.Transport: momapy_bel.core.Reaction,
    momapy.celldesigner.core.Transcription: momapy_bel.core.Reaction,
    momapy.celldesigner.core.Translation: momapy_bel.core.Reaction,
    momapy.celldesigner.core.Modulator: momapy_bel.core.Regulates,
    momapy.celldesigner.core.UnknownModulator: momapy_bel.core.Regulates,
    momapy.celldesigner.core.Inhibitor: momapy_bel.core.Decreases,
    momapy.celldesigner.core.PhysicalStimulator: momapy_bel.core.DirectlyIncreases,
    momapy.celldesigner.core.Catalyzer: momapy_bel.core.DirectlyIncreases,
    momapy.celldesigner.core.Trigger: momapy_bel.core.Increases,
    momapy.celldesigner.core.UnknownCatalyzer: momapy_bel.core.Increases,
    momapy.celldesigner.core.UnknownInhibitor: momapy_bel.core.Decreases,
    momapy.celldesigner.core.Modulation: momapy_bel.core.Regulates,
    momapy.celldesigner.core.Catalysis: momapy_bel.core.DirectlyIncreases,
    momapy.celldesigner.core.Inhibition: momapy_bel.core.Decreases,
    momapy.celldesigner.core.PhysicalStimulation: momapy_bel.core.DirectlyIncreases,
    momapy.celldesigner.core.Triggering: momapy_bel.core.Increases,
    momapy.celldesigner.core.NegativeInfluence: momapy_bel.core.Decreases,
    momapy.celldesigner.core.PositiveInfluence: momapy_bel.core.Increases,
    momapy.celldesigner.core.UnknownModulation: momapy_bel.core.Regulates,
    momapy.celldesigner.core.UnknownCatalysis: momapy_bel.core.DirectlyIncreases,
    momapy.celldesigner.core.UnknownInhibition: momapy_bel.core.Decreases,
    momapy.celldesigner.core.UnknownPositiveInfluence: momapy_bel.core.Increases,
    momapy.celldesigner.core.UnknownNegativeInfluence: momapy_bel.core.Decreases,
    momapy.celldesigner.core.UnknownPhysicalStimulation: momapy_bel.core.DirectlyIncreases,
    momapy.celldesigner.core.UnknownTriggering: momapy_bel.core.Increases,
}

CD_CLASS_TO_PREFERRED_CD_NAMESPACES = {
    momapy.celldesigner.core.Degraded: [],
    momapy.celldesigner.core.SimpleMolecule: SIMPLE_CHEMICAL_NAMESPACES,
    momapy.celldesigner.core.Unknown: [],
    momapy.celldesigner.core.Drug: DRUG_NAMESPACES,
    momapy.celldesigner.core.Ion: SIMPLE_CHEMICAL_NAMESPACES,
    momapy.celldesigner.core.Gene: GENE_NAMESPACES,
    momapy.celldesigner.core.RNA: RNA_NAMESPACES,
    momapy.celldesigner.core.AntisenseRNA: RNA_NAMESPACES,
    momapy.celldesigner.core.Complex: COMPLEX_NAMESPACES,
    momapy.celldesigner.core.GenericProtein: PROTEIN_NAMESPACES,
    momapy.celldesigner.core.Receptor: PROTEIN_NAMESPACES,
    momapy.celldesigner.core.IonChannel: PROTEIN_NAMESPACES,
    momapy.celldesigner.core.TruncatedProtein: PROTEIN_NAMESPACES,
    momapy.celldesigner.core.Phenotype: PHENOTYPE_NAMESPACES,
}


def cd_to_bel(input_file_path, output_file_path):
    cd_map = momapy.io.read(input_file_path)
    cd_model = cd_map.model
    cd_annotations = cd_map.map_element_to_annotations
    bel_model, bel_annotations = cd_model_to_bel_model(
        cd_model, cd_annotations
    )
    momapy.io.write(
        map_=bel_model,
        output_file_path=output_file_path,
        writer="bel",
        annotations=bel_annotations,
    )
    return bel_model, bel_annotations


def cd_model_to_bel_model(cd_model, cd_annotations):
    cd_element_to_bel_element = {}
    bel_model = momapy_bel.core.BELModelBuilder()
    bel_annotations = collections.defaultdict(set)
    for cd_species in cd_model.species:
        _ = _make_and_add_bel_abundance_from_cd_species(
            cd_species,
            cd_annotations,
            bel_model,
            bel_annotations,
            cd_element_to_bel_element,
        )
    for cd_reaction in cd_model.reactions:
        _ = _make_and_add_bel_reaction_from_cd_reaction(
            cd_reaction,
            cd_annotations,
            bel_model,
            bel_annotations,
            cd_element_to_bel_element,
        )
    for cd_modulation in cd_model.modulations:
        _ = _make_and_add_bel_relations_from_cd_modulation(
            cd_modulation,
            cd_annotations,
            bel_model,
            bel_annotations,
            cd_element_to_bel_element,
        )
    bel_model = momapy.builder.object_from_builder(bel_model)
    document_annotation = momapy_bel.core.DocumentAnnotation(
        name="test_name", description="test_description"
    )
    bel_annotations[bel_model].add(document_annotation)
    return bel_model, bel_annotations


def _make_and_add_bel_abundance_from_cd_species(
    cd_species,
    cd_annotations,
    bel_model,
    bel_annotations,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    cd_model_element_cls = type(cd_species)
    bel_model_element_cls = CD_CLASS_TO_BEL_CLASS[cd_model_element_cls]
    bel_abundance = bel_model.new_element(bel_model_element_cls)
    cd_namespaces = CD_CLASS_TO_PREFERRED_CD_NAMESPACES[cd_model_element_cls]
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_species, cd_annotations, cd_namespaces
    )
    if main_iri is not None:
        bel_namespace, bel_identifier = (
            make_bel_namespace_and_identifier_from_cd_iri(main_iri)
        )
    else:
        bel_namespace = CD_NAMESPACE
        bel_identifier = cd_species.name
        bel_identifier = make_normalized_identifier(bel_identifier)
    bel_abundance.namespace = bel_namespace
    bel_abundance.identifier = bel_identifier
    if cd_species.compartment is not None:
        cd_compartment = cd_species.compartment
        _ = _make_and_add_bel_location_from_cd_compartment(
            cd_compartment=cd_compartment,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            bel_annotations=bel_annotations,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_species,
            super_bel_element=bel_abundance,
        )
    if hasattr(cd_species, "modifications") and hasattr(
        bel_abundance, "modifications"
    ):
        for cd_modification in cd_species.modifications:
            if (
                cd_modification.state is not None
                and cd_modification.state.value != "empty"
            ):  # to delete second condition
                _ = _make_and_add_bel_modification_from_cd_modification(
                    cd_modification=cd_modification,
                    cd_annotations=cd_annotations,
                    bel_model=bel_model,
                    bel_annotations=bel_annotations,
                    cd_element_to_bel_element=cd_element_to_bel_element,
                    super_cd_element=cd_species,
                    super_bel_element=bel_abundance,
                )
    if hasattr(cd_species, "structural_states") and hasattr(
        bel_abundance, "modifications"
    ):
        for cd_structural_state in cd_species.structural_states:
            _ = _make_and_add_bel_modification_from_cd_structural_state(
                cd_structural_state=cd_structural_state,
                cd_annotations=cd_annotations,
                bel_model=bel_model,
                bel_annotations=bel_annotations,
                cd_element_to_bel_element=cd_element_to_bel_element,
                super_cd_element=cd_species,
                super_bel_element=bel_abundance,
            )
    bel_subunits = []
    if hasattr(cd_species, "subunits"):
        for cd_subunit in cd_species.subunits:
            bel_subunit = _make_and_add_bel_abundance_from_cd_species(
                cd_species=cd_subunit,
                cd_annotations=cd_annotations,
                bel_model=bel_model,
                bel_annotations=bel_annotations,
                cd_element_to_bel_element=cd_element_to_bel_element,
                super_cd_element=cd_species,
                super_bel_element=bel_abundance,
            )
            bel_abundance.members.add(bel_subunit)
            bel_subunits.append(bel_subunit)
    bel_abundance = momapy.builder.object_from_builder(bel_abundance)
    for bel_subunit in bel_subunits:
        bel_has_component = momapy_bel.core.HasComponent(
            source=bel_abundance, target=bel_subunit
        )
        bel_model.statements.add(bel_has_component)
    if super_cd_element is None:  # not a subunit
        bel_model.statements.add(bel_abundance)
    cd_element_to_bel_element[cd_species] = bel_abundance
    return bel_abundance


def _make_and_add_bel_location_from_cd_compartment(
    cd_compartment,
    cd_annotations,
    bel_model,
    bel_annotations,
    cd_element_to_bel_element,
    super_cd_element,
    super_bel_element,
):
    bel_location = cd_element_to_bel_element.get(cd_compartment)
    if bel_location is None:
        bel_location = bel_model.new_element(momapy_bel.core.Location)
        main_iri = get_main_annotation_iri_from_cd_element(
            cd_compartment,
            cd_annotations,
            cd_namespaces=COMPARTMENT_NAMESPACES,
        )
        if main_iri is not None:
            bel_namespace, bel_identifier = (
                make_bel_namespace_and_identifier_from_cd_iri(main_iri)
            )
        else:
            bel_namespace = CD_NAMESPACE
            bel_identifier = (
                cd_compartment.name if cd_compartment.name else "default"
            )
            bel_identifier = make_normalized_identifier(bel_identifier)
        bel_location.namespace = bel_namespace
        bel_location.identifier = bel_identifier
        cd_element_to_bel_element[cd_compartment] = bel_location
    super_bel_element.location = bel_location
    return bel_location


def _make_and_add_bel_modification_from_cd_modification(
    cd_modification,
    cd_annotations,
    bel_model,
    bel_annotations,
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


def _make_and_add_bel_modification_from_cd_structural_state(
    cd_structural_state,
    cd_annotations,
    bel_model,
    bel_annotations,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_modification = bel_model.new_element(
        momapy_bel.core.ProteinModification
    )
    bel_modification.identifier = make_normalized_identifier(
        cd_structural_state.value
    )
    bel_modification.namespace = CD_NAMESPACE
    bel_modification = momapy.builder.object_from_builder(bel_modification)
    super_bel_element.modifications.append(bel_modification)
    return bel_modification


def _make_and_add_bel_reaction_from_cd_reaction(
    cd_reaction,
    cd_annotations,
    bel_model,
    bel_annotations,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_reaction = bel_model.new_element(momapy_bel.core.Reaction)
    for cd_reactant in cd_reaction.reactants:
        cd_reactant_species = cd_reactant.referred_species
        bel_reactant = cd_element_to_bel_element[cd_reactant_species]
        bel_reaction.reactants.add(bel_reactant)
    for cd_product in cd_reaction.products:
        cd_product_species = cd_product.referred_species
        bel_product = cd_element_to_bel_element[cd_product_species]
        bel_reaction.products.add(bel_product)
    bel_reaction = momapy.builder.object_from_builder(bel_reaction)
    bel_model.statements.add(bel_reaction)
    cd_element_to_bel_element[cd_reaction] = bel_reaction
    for cd_modifier in cd_reaction.modifiers:
        _ = _make_and_add_bel_relations_from_cd_modifier(
            cd_modifier=cd_modifier,
            cd_annotations=cd_annotations,
            bel_model=bel_model,
            bel_annotations=bel_annotations,
            cd_element_to_bel_element=cd_element_to_bel_element,
            super_cd_element=cd_reaction,
            super_bel_element=bel_reaction,
        )
    return bel_reaction


def _make_and_add_bel_relations_from_cd_modifier(
    cd_modifier,
    cd_annotations,
    bel_model,
    bel_annotations,
    cd_element_to_bel_element,
    super_cd_element,
    super_bel_element,
):
    bel_relations = []
    cd_model_element_cls = type(cd_modifier)
    bel_model_element_cls = CD_CLASS_TO_BEL_CLASS[cd_model_element_cls]
    cd_source = cd_modifier.referred_species
    main_iri = get_main_annotation_iri_from_cd_element(
        super_cd_element, cd_annotations, PUBLICATION_NAMESPACES
    )
    if main_iri is not None:
        bel_namespace, bel_identifier = (
            make_bel_namespace_and_identifier_from_cd_iri(main_iri)
        )
        bel_citation = momapy_bel.core.Citation(
            namespace=bel_namespace, identifier=bel_identifier
        )
    else:
        bel_citation = None
    if isinstance(cd_source, momapy.celldesigner.core.Species):
        bel_abundance = cd_element_to_bel_element[cd_source]
        bel_sources = [bel_abundance]
    elif isinstance(cd_source, momapy.celldesigner.core.OrGate):
        bel_sources = []
        for cd_input in cd_source.inputs:
            bel_abundance = cd_element_to_bel_element[cd_input]
            bel_sources.append(bel_abundance)
    elif isinstance(cd_source, momapy.celldesigner.core.AndGate):
        bel_abundance = bel_model.new_element(
            momapy_bel.core.CompositeAbundance
        )
        for cd_input in cd_source.inputs:
            bel_member_abundance = cd_element_to_bel_element[cd_input]
            bel_abundance.members.add(bel_member_abundance)
        bel_abundance = momapy.builder.object_from_builder(bel_abundance)
        bel_sources = [bel_abundance]
    for bel_source in bel_sources:
        bel_relation = bel_model_element_cls(
            source=bel_source, target=super_bel_element
        )
        bel_model.statements.add(bel_relation)
        bel_relations.append(bel_relation)
        if bel_citation is not None:
            bel_annotations[bel_relation].add(bel_citation)
    return bel_relations


def _make_and_add_bel_relations_from_cd_modulation(
    cd_modulation,
    cd_annotations,
    bel_model,
    bel_annotations,
    cd_element_to_bel_element,
    super_cd_element=None,
    super_bel_element=None,
):
    bel_relations = []
    cd_model_element_cls = type(cd_modulation)
    bel_model_element_cls = CD_CLASS_TO_BEL_CLASS[cd_model_element_cls]
    bel_target = cd_element_to_bel_element[cd_modulation.target]
    main_iri = get_main_annotation_iri_from_cd_element(
        cd_modulation, cd_annotations, PUBLICATION_NAMESPACES
    )
    if main_iri is not None:
        bel_namespace, bel_identifier = (
            make_bel_namespace_and_identifier_from_cd_iri(main_iri)
        )
        bel_citation = momapy_bel.core.Citation(
            namespace=bel_namespace, identifier=bel_identifier
        )
    else:
        bel_citation = None
    cd_source = cd_modulation.source
    if isinstance(cd_source, momapy.celldesigner.core.Species):
        bel_abundance = cd_element_to_bel_element[cd_source]
        bel_sources = [bel_abundance]
    elif isinstance(cd_source, momapy.celldesigner.core.OrGate):
        bel_sources = []
        for cd_input in cd_source.inputs:
            bel_abundance = cd_element_to_bel_element[cd_input]
            bel_sources.append(bel_abundance)
    elif isinstance(cd_source, momapy.celldesigner.core.AndGate):
        bel_abundance = bel_model.new_element(
            momapy_bel.core.CompositeAbundance
        )
        for cd_input in cd_source.inputs:
            bel_member_abundance = cd_element_to_bel_element[cd_input]
            bel_abundance.members.add(bel_member_abundance)
        bel_abundance = momapy.builder.object_from_builder(bel_abundance)
        bel_sources = [bel_abundance]
    for bel_source in bel_sources:
        bel_relation = bel_model_element_cls(
            source=bel_source, target=bel_target
        )
        bel_model.statements.add(bel_relation)
        bel_relations.append(bel_relation)
        if bel_citation is not None:
            bel_annotations[bel_relation].add(bel_citation)
    return bel_relations


def get_main_annotation_iri_from_cd_element(
    cd_element, cd_annotations, cd_namespaces
):
    if cd_element in cd_annotations:
        for candidate_cd_namespace in cd_namespaces:
            for cd_annotation in cd_annotations[cd_element]:
                cd_resources = cd_annotation.resources
                for cd_resource in cd_resources:
                    cd_namespace, _ = get_namespace_and_identifier_from_cd_iri(
                        cd_resource
                    )
                    if candidate_cd_namespace == cd_namespace:
                        return cd_resource
    return None


def get_namespace_and_identifier_from_cd_iri(cd_iri):
    cd_iri_split = cd_iri.split(":")
    cd_namespace = cd_iri_split[2]
    cd_identifier = cd_iri_split[3]
    return cd_namespace, cd_identifier


def make_bel_namespace_and_identifier_from_cd_iri(cd_iri):
    cd_namespace, cd_identifier = get_namespace_and_identifier_from_cd_iri(
        cd_iri
    )
    cd_identifier = cd_identifier.replace("%3A", ":")
    get_label_function = CD_NAMESPACE_TO_GET_LABEL_FUNCTION.get(cd_namespace)
    if get_label_function is not None:
        bel_identifier = get_label_function(cd_identifier)
    else:
        bel_identifier = cd_identifier
    if bel_identifier is None:
        bel_identifier = cd_identifier
    bel_identifier = make_normalized_identifier(bel_identifier)
    bel_namespace = CD_NAMESPACE_TO_BEL_NAMESPACE.get(cd_namespace)
    if bel_namespace is None:
        bel_namespace = cd_namespace
    bel_namespace = bel_namespace.upper()
    return bel_namespace, bel_identifier


def make_normalized_identifier(identifier):
    if not identifier.isalnum():
        identifier = f'"{identifier}"'
    return identifier


CD_NAMESPACE_TO_BEL_NAMESPACE = {
    "obo.go": "go",
    "hgnc.symbol": "hgnc",
    "obo.chebi": "chebi",
    "ec-code": "eccode",
    "pubchem.compound": "pubchem",
}

CD_NAMESPACE_TO_GET_LABEL_FUNCTION = {
    "obo.go": cd2bel.utils.get_go_label_from_id,
    "obo.chebi": cd2bel.utils.get_chebi_label_from_id,
    "reactome": cd2bel.utils.get_reactome_label_from_id,
    "interpro": cd2bel.utils.get_interpro_label_from_id,
    "mesh": cd2bel.utils.get_mesh_label_from_id,
    "pfam": cd2bel.utils.get_pfam_label_from_id,
}
