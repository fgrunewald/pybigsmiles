import numpy as np
import networkx as nx
from pysmiles import read_smiles
from .weights import ELEMENTS_TO_WEIGHTS

def compute_weight(graph):
    elements = nx.get_node_attributes(graph, "element")
    weight = 0
    for ele in elements.values():
        weight += ELEMENTS_TO_WEIGHTS[ele]
    return weight

def find_token_indices(line, target):
    idxs = [idx for idx, token in enumerate(line) if token == target]
    for idx in idxs:
        yield idx

def compatible(left, right):
    if left == right:
        return True
    if left[0] == "<" and right[0] == ">":
        if left[1:] == right[1:]:
            return True
    if left[0] == ">" and right[0] == "<":
        if left[1:] == right[1:]:
            return True
    return False

def find_compatible_pair(polymol, residue, bond_type="bond_type", eligible_nodes=None):
    ref_nodes = nx.get_node_attributes(polymol, bond_type)
    target_nodes = nx.get_node_attributes(residue, bond_type)
    for ref_node in ref_nodes:
        if eligible_nodes and\
           polymol.nodes[ref_node]['resid'] not in eligible_nodes:
            continue
        for target_node in target_nodes:
            if compatible(ref_nodes[ref_node],
                          target_nodes[target_node]):
                return ref_node, target_node
    return None

def append_residue(polymol, residue, connection, max_resid):
    mapping = {node: len(polymol.nodes)+node for node in residue.nodes}
    polymol.add_nodes_from([(mapping[node], data) for node, data in residue.nodes(data=True)])
    for node in mapping.values():
        polymol.nodes[node]['resid'] = max_resid + 1
    for edge in residue.edges:
        polymol.add_edge(mapping[edge[0]], mapping[edge[1]])
    if connection:
        for _type in ['bond_type', 'ter_bond_type']:
            if _type in polymol.nodes[connection[0]]:
                del polymol.nodes[connection[0]][_type]
                del polymol.nodes[mapping[connection[1]]][_type]

        polymol.add_edge(connection[0], mapping[connection[1]])
    return polymol

class StochasticObject:
    """
    """
    def __init__(self, residues, end_groups, bonding_probs=None):
        print(residues)
        self.residues = residues
        self.end_groups = end_groups
        self.bonding_probs = bonding_probs
        self.named_objects = {}
        # some useful metrics
        self.resnames = list(self.residues.keys())
        self.end_group_ids = np.arange(0, len(self.end_groups))
        self.mol_weights = {}

        for idx, res in self.residues.items():
            self.mol_weights[idx] = compute_weight(res)
        self.end_weights = {}
        for idx, res in enumerate(self.end_groups):
            self.end_weights[idx] = compute_weight(res)

    def _expand_residue(self, residue):
        mol_graph = residue.copy()
        for node in residue.nodes:
            fragment_name = residue.nodes[node].get('replace', None)
            if fragment_name:
                # we need to generate a completely new random fragment
                if fragment_name in self.named_objects:
                    target_weight = residue.nodes[node].get('target_weight', None)
                    assert target_weight
                    addition = self.named_objects[fragment_name].generate_random(1, target_weight)
                # we simply add an already defined residue fragment
                elif fragment_name in self.residues:
                    addition = self.residues[fragment_name]
                else:
                # no named fragment found; that's not allowed
                    msg = f"Cannot replace fragment with name {fragment_name}."
                    raise ValueError(msg)

                mol_graph.add_nodes_from(addition.nodes(data=True))
                mol_graph.add_edges_from(addition.edges)
                # the smile format defines that the neighboring nodes
                # are uniquely defined and have to connect to the expandable
                # fragment; the fragment itself can have other terminal ends
                # but at this stage they are already capped. Now we assume
                # that the lower index node connects to the frist and the
                # higher index node connects to the second part of the intserted
                # fragment. The actual identifier for terminals is ignored as
                # it cannot be matched against any meanigful surrounding identifier
                # this assumption is due to a ambiguity in the BigSmile definition
                ter_nodes = nx.get_node_attributes(addition, "ter_bond_type")
                assert len(ter_nodes) == 2
                for neighbor, ter in zip(ter_nodes.keys(), residue.neighbors(node)):
                    mol_graph.add_edge(neighbor, ter)

        return mol_graph

    def _add_residue(self, polymol, weight, probabilities, max_resid):
        resname = np.random.choice(self.resnames, probabilities)
        residue = self._expand_residue(self.residues[resname])
        result = find_compatible_pair(polymol, residue)
        if result:
            polymol = append_residue(polymol, residue, result, max_resid)
            weight = weight + self.mol_weights[resname]
        return weight

    def _terminate(self, polymol, weight, max_resid):
        if self.end_groups:
            end_group_id = np.random.choice(self.end_group_ids)
            end_group = self.end_groups[end_group_id]
            result = find_compatible_pair(polymol, end_group, bond_type="ter_bond_type")
            if result:
                polymol = append_residue(polymol, end_group, result, max_resid)
                weight = weight + self.end_weights[end_group_id]
        return polymol, weight

    def generate_random(self,
                        max_resid=0,
                        random_terminate=False,
                        target_weight=None):
        """
        """
        polymol = nx.Graph()
        probabilities = None
        resid = np.random.choice(list(self.residues.keys()),
                                 probabilities)
        residue = self.residues[resid]
        polymol.add_nodes_from(residue.nodes(data=True))
        polymol.add_edges_from(residue.edges)
        weight = self.mol_weights[resid]
        polymol, weight = self._terminate(polymol, weight, max_resid)
        print(polymol.nodes(data=True))
        while True:
            weight = self._add_residue(polymol,
                                       weight,
                                       probabilities,
                                       max_resid)
            if random_terminate and np.random.uniform(0, 1) > 0.5:
                polymol, weight = self._terminate(polymol, weight, max_resid)
                break
            if weight >= target_weight:
                break
        polymol, weight = self._terminate(polymol, weight, max_resid)
        return polymol

    def generate_molecules(self, n, target_weight):
        for idx in range(0, n):
            yield self.generate_random(target_weight=target_weight)

    @classmethod
    def from_bigsmile(cls, line):
        """
        Create a stochastic object from a big smile.
        Tdo:
        Raise error when endgroup is missing []
        """
        line = line.strip()
        # a stochastic object is enclosed in '{' and '}'
        start_idx = next(find_token_indices(line, "{"))
        stop_idx = next(find_token_indices(line, "}"))
        stoch_line = line[start_idx+1:stop_idx]
        # residues are separated by , and end
        # groups by ;
        if ';' in stoch_line:
            residue_string, terminii_string = stoch_line.split(';')
        else:
            residue_string = stoch_line
            terminii_string = None
        # let's read the smile residue strings
        residues = []
        count = 0
        for residue_string in residue_string.split(','):
            # figure out if this is a named object
            if residue_string[0] == "#":
                jdx = next(find_token_indices(residue_string, "="))
                name = residue_string[:jdx]
                residue_string = residue_string[jdx:]
            else:
                name = count

            mol_graph = read_smiles(residue_string)
            residues.append((name, mol_graph))
            count += 1
        # let's read the terminal residue strings
        end_groups = []
        if terminii_string:
            for terminus_string in terminii_string.split(','):
                mol_graph = read_smiles(terminus_string)
                bond_types = nx.get_node_attributes(mol_graph, "bond_type")
                nx.set_node_attributes(mol_graph, bond_types, "ter_bond_type")
                end_groups.append(mol_graph)
        return cls(dict(residues), end_groups)
