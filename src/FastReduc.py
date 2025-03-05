# Copyright 2025 OTH - Laboratory for Digitalisation (LfD)
# Written by Lukas Schmidbauer
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

import pprint
import networkx as nx
from itertools import combinations
from sortedcontainers import SortedDict


def fastPolyQuadratisation(
    poly_dict,
    max_degree,
    reduction_variable_prefix=None,
    var_pair_choice=None,
    selection_quantile=1.0,
):
    """
    Quadratises a binary polynomial with efficient data structures in the background.
    Works on compact polynomials (i.e. {(1,): 2, (0,1): 3)}. Only numeric variable names are allowed.

    Input:  poly_dict: a quark PolyBinary dictionary
            max_degree: maximum degree of the resulting polynomial (typically 2)
            reduction_variable_prefix: Prefix for new variables (Not implemented yet)
            var_pair_choice: Strategy for selecting the next variable pair (Not implemented -> use selection_quantile)
            selection_quantile: value between 0 and 1. 0: sparse, 1: dense

    Output: Reduced PolyBinary and penalty terms to be able to scale them separately
    """

    monomial_idx_dict = dict()
    i = 0
    for monomial in poly_dict.keys():
        monomial_idx_dict[i] = monomial
        i += 1
    
    G = nx.MultiGraph()
    G.add_nodes_from(poly_dict.variables)

    for idx in monomial_idx_dict.keys():
        weightless_edges = list(combinations(monomial_idx_dict[idx], 2))
        for ed in weightless_edges:
            G.add_edge(ed[0], ed[1], key=idx)

    sorted_ranking_dict = SortedDict()
    for node in G.nodes():
        for neighbour in G.neighbors(node):
            cnt = G.number_of_edges(node, neighbour)
            if cnt > 1 and node < neighbour:
                if sorted_ranking_dict.__contains__(cnt):
                    sorted_ranking_dict[cnt][(node, neighbour)] = ""
                else:
                    sorted_ranking_dict[cnt] = dict()
                    sorted_ranking_dict[cnt][(node, neighbour)] = ""

    poly_dict = dict(poly_dict)

    penalty_dict = dict()

    current_degree = max(map(len, poly_dict.keys()), default=0)

    NewVarIdx = 0
    for m in poly_dict.keys():
        for idx in m:
            if (
                idx > NewVarIdx
            ): 
                NewVarIdx = idx
    NewVarIdx += 1

    while current_degree > max_degree:
        if sorted_ranking_dict.__len__() == 0:
            break

        sdidx = int((sorted_ranking_dict.__len__() - 1) * selection_quantile)
        var_pair_choice = next(iter(sorted_ranking_dict.peekitem(sdidx)[1]))
        

        xi = var_pair_choice[0]
        xj = var_pair_choice[1]

        c4_monomials = set()
        for k in G[xi][xj].keys():
            c4_monomials.add(k)

        c4_edges_old = set()
        c4_edges_new = set()
        for monomial_idx in c4_monomials:
            monomial_old = monomial_idx_dict[monomial_idx]
            weightless_edges = list()
            weightless_edges.append([xi, xj])
            for var in monomial_old:
                if (var != xi) and (var != xj):
                    weightless_edges.append([var, xi])
                    weightless_edges.append([var, xj]) 
            for ed in weightless_edges:
                if ed[0] < ed[1]:
                    c4_edges_old.add((ed[0], ed[1], monomial_idx))
                else:
                    c4_edges_old.add((ed[1], ed[0], monomial_idx))

            monomial_new = replaceVariablePair((xi,xj), NewVarIdx, monomial_old)
            weightless_edges = list()
            for var in monomial_new:
                if var != NewVarIdx:
                    weightless_edges.append([var, NewVarIdx])
            for ed in weightless_edges:
                if ed[0] < ed[1]:
                    c4_edges_new.add((ed[0], ed[1], monomial_idx))
                else:
                    c4_edges_new.add((ed[1], ed[0], monomial_idx))

        for monomial in c4_monomials:
            m = monomial_idx_dict[monomial]
            mnew = replaceVariablePair(var_pair_choice, NewVarIdx, m)
            monomial_idx_dict[monomial] = mnew
            poly_dict[mnew] = poly_dict[m]
            poly_dict.__delitem__(m)
        penalty_dict[(NewVarIdx,)] = 3
        penalty_dict[(var_pair_choice[0], var_pair_choice[1])] = 1
        penalty_dict[(NewVarIdx, var_pair_choice[0])] = -2
        penalty_dict[(NewVarIdx, var_pair_choice[1])] = -2

        
        P_lost_node_pairs = set()
        for edge in c4_edges_old:
            if edge[0] < edge[1]:
                P_lost_node_pairs.add((edge[0], edge[1], G.number_of_edges(edge[0], edge[1])))
            else:
                P_lost_node_pairs.add((edge[1], edge[0], G.number_of_edges(edge[0], edge[1])))
        
        P_new_node_pairs = set()
        for edge in c4_edges_new:
            if edge[0] < edge[1]:
                P_new_node_pairs.add((edge[0], edge[1], G.number_of_edges(edge[0], edge[1])))
            else:
                P_new_node_pairs.add((edge[1], edge[0], G.number_of_edges(edge[0], edge[1])))

        for edge in c4_edges_old:
            G.remove_edge(edge[0], edge[1], key=edge[2])
        for edge in c4_edges_new:
            G.add_edge(edge[0], edge[1], key=edge[2])

        for np in P_lost_node_pairs.union(P_new_node_pairs):
            if np[0] > np[1]:
                np = (np[1], np[0], np[2])
            if np[2] > 1:
                sorted_ranking_dict[np[2]].__delitem__((np[0], np[1]))
            cnt = G.number_of_edges(np[0], np[1])
            if cnt > 1:
                if sorted_ranking_dict.__contains__(cnt):
                    sorted_ranking_dict[cnt][(np[0], np[1])] = ""
                else:
                    sorted_ranking_dict[cnt] = dict()
                    sorted_ranking_dict[cnt][(np[0], np[1])] = ""

            if sorted_ranking_dict.__contains__(np[2]):
                if sorted_ranking_dict[np[2]].__len__() == 0:
                    sorted_ranking_dict.__delitem__(np[2])

        NewVarIdx += 1

    current_degree = max(map(len, poly_dict.keys()), default=0)

    if current_degree > max_degree:
        unreduced_monomials = list(poly_dict.keys())

        for m in unreduced_monomials:
            m_degree = len(m)
            mnew = []
            mlastit = m
            while m_degree > 2:
                for i in range(0, m_degree - 1, 2):
                    choice = (mlastit[i], mlastit[i + 1])
                    mnew.append(NewVarIdx)

                    penalty_dict[(NewVarIdx,)] = 3
                    penalty_dict[(choice[0], choice[1])] = 1
                    penalty_dict[(NewVarIdx, choice[0])] = -2
                    penalty_dict[(NewVarIdx, choice[1])] = -2

                    NewVarIdx += 1
                if m_degree % 2 == 1:
                    mnew.append(mlastit[-1])

                m_degree = len(mnew)
                mlastit = mnew
                mnew = []

            alpha = poly_dict[m]
            poly_dict.__delitem__(m)
            poly_dict[tuple(mlastit)] = alpha

    return poly_dict, penalty_dict

def replaceVariablePair(pair, newVarIdx, monomial):
    """
    Replaces the variable pair in monomial and returns the new result
    Input:
            pair: variable pair (e.g., (1, 2))
            newVarIdx: idx of the new variable (e.g. 71787)
            monomial: concrete monomial (e.g., (1,2,4,6))
    """
    mnew = []
    for i in range(len(monomial)):
        if pair[0] == monomial[i] or pair[1] == monomial[i]:
            pass
        else:
            mnew.append(monomial[i])
    mnew.append(newVarIdx)
    return tuple(mnew)
