import enum
import networkx as nx
from .read_smiles import read_smiles
from .stocastic_object import StochasticObject, append_residue, find_compatible_pair
from .gsmile_helpers import _parse_non_st_obj_gsmile

@enum.unique
class BigTokenType(enum.Enum):
    """Possible BigSMILES token types"""
    START = 1
    STOP = 2
    APPEND = 3
    SMILE = 4
    GSMILE_OPEN = 5
    GSMILE_CLOSE = 6
    GSMILE_APPEND = 7

class _Smile:

    def __init__(self, smile_str):
        self.graph = read_smiles(smile_str)

    def generate_random(self,
                        target_weight=None):
        return self.graph

    def generate_seq(self, sequence=None):
        return self.graph

class BigPolymer(nx.DiGraph):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.preprocessed = False
        self.max_resid = 1
        self.gsmile_attributes = {}

    def preprocess(self):
        if self.preprocessed:
            print("BigPolymer is already preprocessed. Skipping it.")
        else:
            for node in self.nodes:
                if "obj_str" in self.nodes[node]:
                    obj_str = self.nodes[node]["obj_str"]
                    st_obj = StochasticObject.from_bigsmile(obj_str)
                    self.nodes[node]["st_obj"] = st_obj
                elif "smile" in self.nodes[node]:
                    self.nodes[node]["residue"] = _Smile(self.nodes[node]["smile"])
                else:
                    raise IOError
            self.preprocessed = True

    def generate_random(self):
        if not self.preprocessed:
            self.preprocess()

        polymol = nx.Graph()
        max_resid = 0
        residue = self.nodes[0]["residue"].generate_random(self.nodes[0]["weight"])
        polymol = append_residue(polymol, residue, None, max_resid)
        for prev_node, current_node in nx.bfs_edges(self, source=0):
            residue = self.nodes[current_node]["residue"].generate_random(self.nodes[current_node]["weight"])
            result = find_compatible_pair(polymol, residue, eligible_nodes=prev_node_resids)
            if result:
                polymol = append_residue(polymol, residue, result, max_resid)
                max_resid += 1
            else:
                # there has to be at least one connectable node
                # even of mixed terminal definitions are used
                raise IOError

class BigSmileParser:

    def __init__(self):
        self.current_object = None
        self.current_smile = ""
        self.current_gsmile = ""
        self.unfinished_objects = []
        self.smiles = []
        self.node_count = 0
        self.ref_node = None
        self.nest_counter = 0
        self.parser_actions = {BigTokenType.START: self.new_stochastic_object,
                               BigTokenType.STOP: self.finish_stochastic_object,
                               BigTokenType.APPEND: self.append_stochastic_object,
                               BigTokenType.SMILE: self.append_smile,
                               BigTokenType.GSMILE_OPEN: self.open_gsmile,
                               BigTokenType.GSMILE_CLOSE: self.close_gsmile,
                               BigTokenType.GSMILE_APPEND: self.append_gsmile}

        self.polymer = BigPolymer()

    def new_stochastic_object(self, token):
        if self.current_smile:
            self.polymer.add_node(self.node_count,
                                  smile_str = self.current_smile)
            self.polymer.add_edge(self.ref_node, self.node_count)
            self.ref_node = self.node_count
            self.node_count += 1
            self.current_smile = ""

        if self.current_object is not None:
            # here we need to take care of the nesting
            self.current_object += f"#{nest_count}"
            self.unfinished_objects.append(self.current_object)
        self.current_object = f"#{nest_count}={token}"

    def finish_stochastic_object(self, token):
        self.current_object += token
        # create a new node and store reference
        self.polymer.add_node(self.node_count,
                              obj_str=self.current_object)
        self.polymer.add_edge(self.ref_node, self.node_count)
        self.ref_node = self.node_count
        self.node_count += 1
        if self.unfinished_objects:
            self.current_object = self.unfinished_objects[-1]
        else:
            self.current_object = None

    def open_gsmile(self, token):
        if self.current_smile:
            self.polymer.add_node(self.node_count,
                                  smile_str = self.current_smile)
            self.node_count += 1
            self.current_smile = ""
        self.current_gsmile = token

    def close_gsmile(self, token):
        key, value = _parse_non_st_obj_gsmile(self.current_gsmile)
        self.polymer.gsmile_attributes[key] = value
        self.current_gsmile = ""

    def append_stochastic_object(self, token):
        self.current_object += token

    def append_smile(self, token):
        self.current_smile += token

    def append_gsmile(self, token):
        self.current_gsmile += token

    def tokenize(self, line):
        line_iter = iter(line)
        for token in line_iter:
            # start of a stochastic object
            if token == "{":
                token_type = BigTokenType.START
            # end of a stochastic object
            elif token == "}":
                token_type = BigTokenType.STOP
            # we have an open stochastic object and
            # append
            elif self.current_object:
                token_type = BigTokenType.APPEND
            # here we finish parsing a generative
            # bigsmile string that is not included
            # in the stochastic object
            elif token == "|" and self.current_gsmile:
                token_type = BigTokenType.GSMILE_CLOSE
            # here we start parsing a generative
            # big smile definition
            elif token == "|" and not self.current_gsmile:
                token_type = BigTokenType.GSMILE_CLOSE
            # here we append a generative big smile
            # definition
            elif self.current_gsmile:
                token_type = BigTokenType.GSMILE_APPEND
            # if it is not a stochastic object and it is
            # not a generative bigsmile string it must be
            # a regular end-group smile string
            else:
                token_type = BigTokenType.SMILE

            yield token, token_type

    def finish(self):
        if self.current_smile:
            self.polymer.add_node(self.node_count,
                                  smile_str=self.current_smile)
            self.polymer.add_edge(self.ref_node, self.node_count)
            self.ref_node = self.node_count
            self.node_count += 1
            self.current_smile = ""

    def parse(self, big_smile_line):
        """
        Parse big smile objects.
        """
        for token, token_type in self.tokenize(big_smile_line):
            self.parser_actions[token_type](token)
        self.finish(self)
