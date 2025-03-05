# Copyright 2025 <>
# Written by <>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import random
import sympy as sp
import numpy as np
import math
import scipy.special
import networkx as nx
from multiset import Multiset
from copy import deepcopy
from quark import PolyBinary
from itertools import combinations


from FastReduc import *

def createSATinstance(k, num_clauses, num_variables, seed = 42):
    """
    Creates a random exact k-SAT instance from the set [1;num_variables] variables and with exactly [1;num_clauses] clauses randomly.
    Note that not all variables might be used in literals
    """
    random.seed(seed)
    literals = ["x" + str(i + 1) for i in range(num_variables)] + ["!x" + str(i + 1) for i in range(num_variables)]

    formula = []

    for i in range(num_clauses):
        C = random.sample(literals, k)
        formula.append(C)

    return formula

################################################################################################################################################################
# Transformations
################################################################################################################################################################

def ChoisReduction(SAT_instance):
    """
    Creates a QUBO for given SAT_instance by firstly mapping it to MIS
    Input: 
        SAT_instance:   List of lists. Inner lists represent clauses and contain string literals. 
                        A negated literal is represented as "!<variable>"
    """
    G = nx.Graph()
    for num,C in enumerate(SAT_instance):
        nodes = [(num, lit) for lit in C]
        edges = combinations(nodes,2)
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)
    
    for n1 in G.nodes:
        for n2 in G.nodes:
            if (n1[0] != n2[0]) and (n1[1] != n2[1]):
                if n1[1].replace("!", "") == n2[1] or n1[1] == n2[1].replace("!", ""):
                    G.add_edge(n1, n2)

    relable_dict = dict()
    for n in G.nodes:
        relable_dict[n] = str(n[0]) + "c" + str(n[1])
    G = nx.relabel_nodes(G, relable_dict)

    vars = []
    idx = 0
    qubo = np.identity(len(G.nodes))
    for n in G.nodes:
        vars.append(n)

    for r in range(len(vars)):
        for c in range(r, len(vars)):
            if G.has_edge(vars[r],vars[c]):
                qubo[r][c] = 2
                qubo[c][r] = 2

    return G, qubo

class Binary(sp.Symbol):
    def _eval_power(self, other):
        return self


def recursive_expand(input:list, outDict:dict, path=[]):
    """
    Rekursively expands a term of the form (1-x1)(1-x2)...(1-xn) and writes result in outDict in the form {<monomial>: <alpha>}
    Input: List of variables to expand, for example ["x1", "x3", ...]. They are expanded, assuming above form: (1-x1)...
    Care when using this method with PolyBinary! -> Exponents are ignored and terms are added and subtracted there.
    """
    if len(input) <= 0: 
        if len(path) % 2 == 0:
            outDict.update({tuple(path): 1})
        else:
            outDict.update({tuple(path): -1})
    else:
        nv = input[0]
        p1 = path + [nv]
        p2 = path

        recursive_expand(input[1:], outDict, p1)
        recursive_expand(input[1:], outDict, p2)



def PUBOfromSAT(SAT_instance):
    """
        Creates a PUBO from a SAT_instance, via extended DeMorgan's law
    """
    vars = set()
    for C in SAT_instance:
        for lit in C:
            vars.add(lit.replace("!", ""))

    pubo = PolyBinary()
    for C in SAT_instance:
        negativeVars = set()
        positiveVars = set()
        for lit in C:
            if "!" in lit:
                negativeVars.add(int(lit.replace("!x","")))
            else:
                positiveVars.add(int(lit.replace("x", "")))
            
        positiveExpandDict = dict()
        recursive_expand(list(positiveVars), positiveExpandDict)

        puboDict = dict()
        for key in positiveExpandDict.keys():
            puboDict[tuple(sorted(tuple(negativeVars) + key))] = positiveExpandDict[key]
        

        pubo += 1 - PolyBinary(puboDict)
    
    return PolyBinary({}) - pubo


def textbook3SatfromkSat(SAT_instance):
    """
    Iteratively reduces the size of exact k-SAT clauses to exact 3-SAT clauses by introducing new variables:
    For example: (x1 + x2 + x3 + !x4) -> (x1 + x2 + x5) * (!x5 + x3 + !x4)
    (x1 + x2 + x3 + x4 + x5 + x6) -> (x1 + x2 + x3 + x4 + x7) * (!x7 + x5 + x6) -> (x1 + x2 + x3 + x8) * (!x8 + x4 + x7) * (!x7 + x5 + x6) 
        -> (x1 + x2 + x9) * (!x9 + x3 + x8) * (!x8 + x4 + x7) * (!x7 + x5 + x6)
    """
    n = 0
    for C in SAT_instance:
        for lit in C:
            n = max(int(lit.replace("!", "").replace("x", "")), n)

    n += 1

    nSAT = []
    for C in SAT_instance:
        if len(C) <= 3:
            nSAT.append(C)
        else:
            Ctmp = deepcopy(C)
            while len(Ctmp) > 3:
                xi = Ctmp.pop()
                xj = Ctmp.pop()
                Ctmp.append("x" + str(n))
                
                nSAT.append(["!x" + str(n), xj, xi])
                n+=1      
            nSAT.append(Ctmp)

    return nSAT

def QuboMatToPolyDict(QuboMatrix):
    """
        Returns the dictionary representation of a qubo matrix, that is {<monomial>: <alpha>}. For example {(1,2,3): 13} for 13*x1x2x3
    """
    poly = {}
    for d in range(len(QuboMatrix)):
        if QuboMatrix[d][d] != 0:
            poly[(d,)] = QuboMatrix[d][d]

    for r in range(len(QuboMatrix)):
        for c in range(r+1, len(QuboMatrix[r])):
            if (QuboMatrix[r][c] + QuboMatrix[c][r]) != 0:
                poly[(r,c)] = QuboMatrix[r][c] + QuboMatrix[c][r]

    return poly

def PolyDictToQubo(polynomial):
    """
        Returns upper triangular qubo representation of quadratic polynomial (PolyBinary)
    """
    polynomial = polynomial.compact()

    qubo = np.zeros(shape=(len(polynomial.variables), len(polynomial.variables)))
    for m in polynomial.keys():
        if len(m) > 2:
            print("Supplied higher order monomial! -> Ignoring")
        elif len(m) == 2:
            qubo[m[0]][m[1]] = polynomial[m]
        elif len(m) == 1:
            qubo[m[0]][m[0]] = polynomial[m]
        else:
            print("Supplied constant monomial! -> Ignoring")
    
    return qubo

################################################################################################################################################################
#Properties
################################################################################################################################################################

def getSATVars(SAT_instance):
    """
        Returns set of used variables from SAT instance
    """
    vars = set()
    for C in SAT_instance:
        for lit in C:
            vars.add(lit.replace("!", ""))
    return vars

def getSATliterals(SAT_instance):
    """
        Returns set of used literals from SAT_instance
    """
    lits = set()
    for C in SAT_instance:
        for lit in C:
            lits.add(lit)
    return lits

def getDensities(poly):
    """
    Input:
        poly: a dictionary of the form: {(<tuple of variables>): <alpha>, ...}
    Returns all degree-k densities of poly (i.e. actual / possible monomials of degree-k) in the form of a dictionary
        {<k>: <deg-k-density>}
    """

    variables = set()
    degree_dict = dict()
    for monomial in poly.keys():
        for v in monomial:
            variables.add(v)
        current_degree = len(monomial)
        if degree_dict.__contains__(current_degree):
            degree_dict[current_degree].append((poly[monomial], monomial))
        else:
            degree_dict[current_degree] = [(poly[monomial], monomial)]

    densities_dict = dict()
    for i in range(max(len(degree_dict.keys()) + 1, 5)):
        if degree_dict.__contains__(i):
            densities_dict[i] = len(degree_dict[i]) / scipy.special.comb(
                len(variables), i, exact=True
            )
        else:
            densities_dict[i] = 0

    return densities_dict

################################################################################################################################################################
# Visuals
################################################################################################################################################################

def MultiGraphToLatexTikz(G, radius, clique_sort=False):
        """
        Input: networkx MultiGraph with compact node names
                radius: Node radius (typ 2, big graphs 10)
                clique_sort: Sorts nodes by clique size
        Returns: latex string of that drawn MultiGraph
        """

        oG = deepcopy(G)
        nG = nx.Graph()
        nG.add_nodes_from(list(G.nodes()))
        nG.add_edges_from(list(G.edges()))
        
        out = "\\resizebox{\\linewidth}{!}{%\n\\begin{tikzpicture}\n\t\\begin{scope}[circle]\n\t\t\draw\n"
        
        size = G.number_of_nodes()
        i = 0

        node_colour = nx.get_node_attributes(oG, "colour", default="lfdblack")

        if clique_sort:
            nodes = []
            while len(G.nodes) > 0:
                n = sorted(list(nx.find_cliques(G)),key=len, reverse=True)[0]
                nodes += n
                G.remove_nodes_from(n)
                
            for n in nodes:
                out += "\t\t(" + str((i*360.0)/size) +  ":" + str(radius) + ") \tnode [" + node_colour[n] + "] (x" + str(n) + ") {$x_{" + str(n) + " }$}\n"
                i+=1
        else:
            for n in G.nodes():
                out += "\t\t(" + str((i*360.0)/size) +  ":" + str(radius) + ") \tnode [" + node_colour[n] + "] (x" + str(n) + ") {$x_{" + str(n) + " }$}\n"
                i+=1
        
        out += "\t\t;\n\t\\end{scope}\n"
        
        edge_colour = nx.get_edge_attributes(oG, "colour", default="lfdblack")
        edge_weight = nx.get_edge_attributes(oG, "weight", default=1)
        
        out += "\t\\begin{scope}[-]\n"
        for e in nG.edges():
            num_ed = edge_weight[e]
            if num_ed == 1:
                out += "\t\t\\draw [" + edge_colour[e] + "] (x" + str(e[0]) + ") to " + "(x" + str(e[1]) + ");\n"
            else:
                out += "\t\t\\draw [" + edge_colour[e] + "] (x" + str(e[0]) + ") to node[above, sloped] {$" + str(num_ed) + "$} (x" + str(e[1]) + ");\n"
            
        
        out += "\t\\end{scope}\n\end{tikzpicture}\n}"
        
        return out

def genCouplingGraphFromPolynomial(polynomial, variableIdxAssignmentDict):
    """
        Generates a multi graph from a given quark polynomial.
        Variables = nodes, variables in the same monomial = edge between them
        variableIdxAssignement: Dictionary with original variable names as keys (as they appear in polynomial and indices [0,...] as values
    """
    G = nx.MultiGraph()
    nodes = []
    for var in polynomial.variables:
        nodes += [variableIdxAssignmentDict[var]]
    
    G.add_nodes_from(nodes)

    edges = []
    for monomial in polynomial.keys():
        newMonomial = tuple()
        for i in monomial:
            newMonomial += (variableIdxAssignmentDict[i],)
        edges += (list(combinations(newMonomial, 2)))

    G.add_edges_from(edges)

    return G

def primalGraphFromPolynomial(polynomial, oldGraph=None):
    """
        Creates a primal (or variable-incidence) graph from a given polynomial (in Quark format).
        Variables are nodes and edges are introduced whenever two variables occur in the same monomial
        If oldGraph is provided, new edges and new nodes will be coloured red 
    """
    edges = Multiset()
    vars = set()
    for m in polynomial.keys():
        sub_vars = set()
        for var in m:
            vars.add(var)
        for e in list(combinations(m, 2)):
            edges.add(e)
    
    G = nx.Graph()
    for n in vars:
        if oldGraph is not None and n not in oldGraph.nodes:
            G.add_node(n, colour="lfdred")
        else:
            G.add_node(n, colour="lfdblack")

    for e in edges:
        if oldGraph is not None and e not in oldGraph.edges:
            G.add_edge(e[0], e[1], weight=edges.get(e, 0), colour="lfdred")
        else:
            G.add_edge(e[0], e[1], weight=edges.get(e, 0), colour="lfdblack")

    
    return G

def primalGraphFromkSAT(SAT_instance, oldGraph=None):
    """
        Creates a primal (or variable-incidence) graph from a given SAT instance.
        Variables are nodes and edges are introduced whenever two variables occur in the same clause (doesn't matter if negated or not)
        If oldGraph is provided, new edges and new nodes will be coloured red 
    """

    edges = []
    vars = set()
    for C in SAT_instance:
        sub_vars = set()
        for lit in C:
            vars.add(int(lit.replace("!", "").replace("x","")))
            sub_vars.add(int(lit.replace("!", "").replace("x","")))
        edges += combinations(sub_vars, 2)
    
    G = nx.Graph()
    for n in vars:
        if oldGraph is not None and n not in oldGraph.nodes:
            G.add_node(n, colour="lfdred")
        else:
            G.add_node(n, colour="lfdblack")

    for e in edges:
        if oldGraph is not None and e not in oldGraph.edges:
            G.add_edge(e[0], e[1], weight=edges.count(e), colour="lfdred")
        else:
            G.add_edge(e[0], e[1], weight=edges.count(e), colour="lfdblack")

    return G