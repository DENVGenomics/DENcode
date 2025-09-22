import networkx as nx
import json

def export_graph(P_matrix, output_path="output/transmission_graph.graphml"):
    G = nx.DiGraph()
    for src, targets in P_matrix.items():
        for tgt, prob in targets.items():
            if prob > 0:
                G.add_edge(src, tgt, weight=prob)
    nx.write_graphml(G, output_path)

def save_pij_json(P_matrix, output_path="output/pij_matrix.json"):
    with open(output_path, "w") as f:
        json.dump(P_matrix, f, indent=2)
