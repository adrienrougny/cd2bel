# adapted from pybel's source code
def to_neo4j(graph, neo_connection, collection_name=None):
    """Upload a BEL graph to a Neo4j graph database using :mod:`py2neo`.

    :param pybel.BELGraph graph: A BEL Graph
    :param neo_connection: A :mod:`py2neo` connection object. Refer to the
     `py2neo documentation <http://py2neo.org/v3/database.html#the-graph>`_ for how to build this object.
    :type neo_connection: str or py2neo.Graph

    Example Usage:

    >>> import py2neo
    >>> import pybel
    >>> from pybel.examples import sialic_acid_graph
    >>> neo_graph = py2neo.Graph("http://localhost:7474/db/data/")  # use your own connection settings
    >>> pybel.to_neo4j(sialic_acid_graph, neo_graph)
    """
    import py2neo
    from pybel.constants import (
        ANNOTATIONS,
        CITATION,
        EVIDENCE,
        FUSION,
        MEMBERS,
        NAMESPACE,
        RELATION,
        SOURCE_MODIFIER,
        TARGET_MODIFIER,
        VARIANTS,
    )
    from pybel.utils import flatten_dict

    if isinstance(neo_connection, str):
        neo_connection = py2neo.Graph(neo_connection)

    tx = neo_connection.begin()

    document_node = py2neo.Node("Document", name=graph.name)
    tx.create(document_node)

    if collection_name is not None:
        cursor = tx.run(
            f"MERGE (collection:Collection {{name: '{collection_name}'}}) RETURN collection"
        )
        for record in cursor:
            for collection_node in record:
                break
        rel = py2neo.Relationship(
            collection_node, "HAS_DOCUMENT", document_node
        )
        tx.create(rel)

    node_map = {}

    nodes = list(graph)
    for node in nodes:
        if (
            NAMESPACE not in node
            or VARIANTS in node
            or MEMBERS in node
            or FUSION in node
        ):
            attrs = {"name": node.as_bel()}
        else:
            attrs = {"namespace": node.namespace}

            if node.name and node.identifier:
                attrs["name"] = node.name
                attrs["identifier"] = node.identifier
            elif node.identifier and not node.name:
                attrs["name"] = node.identifier
            elif node.name and not node.identifier:
                attrs["name"] = node.name

        py2neo_node = py2neo.Node(node.function, **attrs)
        node_map[node] = py2neo_node

        tx.create(py2neo_node)

        rel = py2neo.Relationship(document_node, "HAS_NODE", py2neo_node)
        tx.create(rel)

    edges = graph.edges(keys=True, data=True)

    for u, v, key, node in edges:
        rel_type = node[RELATION]

        d = node.copy()
        del d[RELATION]

        attrs = {}

        annotations = d.pop(ANNOTATIONS, None)
        if annotations:
            for annotation, values in annotations.items():
                # following code modified by A. Rougny
                non_dict_values = []
                for value in values:
                    if isinstance(value, dict):
                        if "namespace" in value and "identifier" in value:
                            non_dict_values.append(
                                f"{value['namespace']}_{value['identifier']}"
                            )
                        else:
                            print(annotation, value)
                    else:
                        non_dict_values.append(value)
                attrs[annotation] = non_dict_values
                # end of modified code by A. Rougny
        citation = d.pop(CITATION, None)
        if citation:
            attrs[CITATION] = citation.curie

        if EVIDENCE in d:
            attrs[EVIDENCE] = d[EVIDENCE]

        for side in (SOURCE_MODIFIER, TARGET_MODIFIER):
            side_data = d.get(side)
            if side_data:
                attrs.update(flatten_dict(side_data, parent_key=side))

        rel = py2neo.Relationship(
            node_map[u], rel_type, node_map[v], key=key, **attrs
        )
        tx.create(rel)

    tx.commit()
