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

import time
import pandas as pd
import glob
from util import *
from p_tqdm import p_imap

def do_benchmark():
    param_k = []
    param_vars = []
    param_clauses = []
    for k in [3,5,7,9,10,15,20]:
        for vars in range(k, 60, 3):
            for clauses in range(2, 300, 13):
                param_k.insert(0,k)
                param_vars.insert(0,vars)
                param_clauses.insert(0,clauses)

    ########################################################
    # k-Sat -Textbook-> 3-SAT -Choi(MIS)-> QUBO
    ########################################################

    completeKT3CMCQ = p_imap(singleKT3CMCQ, param_k, param_clauses, param_vars)
    pdresKT3CMCQ = pd.DataFrame(completeKT3CMCQ)
    pdresKT3CMCQ.to_csv("Experiments/KT3CMCQ.csv")

    #######################################################
    # k-Sat -Choi*(MIS)-> QUBO
    ########################################################

    completeKCMCQ = p_imap(singleKCMCQ, param_k, param_clauses, param_vars)
    pdresKCMCQ = pd.DataFrame(completeKCMCQ)
    pdresKCMCQ.to_csv("Experiments/KCMCQ.csv")

    ########################################################
    # k-Sat -DeMorgan-> PUBO -Fast-> QUBO
    ########################################################
    
    #extra parameter: percentile for fast quadratisation
    p_percentile = []
    p_k = []
    p_clauses = []
    p_vars = []
    for p in [0.01, 0.05, 0.1, 0.4, 0.7, 1.0]:
        for i in range(len(param_k)):
            p_k.append(param_k[i])
            p_clauses.append(param_clauses[i])
            p_vars.append(param_vars[i])
            p_percentile.append(p)


    completeKDePFQ = p_imap(singleKDePFQ, p_k, p_clauses, p_vars, p_percentile)
    pdresKDePFQ = pd.DataFrame(completeKDePFQ)
    pdresKDePFQ.to_csv("Experiments/KDePFQ.csv")

    ########################################################
    # k-Sat -Textbook-> 3-SAT -Dobrynin-> PUBO -Fast-> QUBO
    ########################################################

    completeKT3DPFQ = p_imap(singleKT3DPFQ, p_k, p_clauses, p_vars, p_percentile)
    pdresKT3DPFQ = pd.DataFrame(completeKT3DPFQ)
    pdresKT3DPFQ.to_csv("Experiments/KT3DPFQ.csv")


def do_graph_eval(clique_sort=False):
    k = 6
    clauses = 5
    vars = 20
    percentile = 1.0

    folder = "GraphEval"

    print("KT3CMCQ")
    singleKT3CMCQ(k,    clauses, vars,              folder)
    print("KCMCQ")
    singleKCMCQ(k,      clauses, vars,              folder)
    print("KDePFQ")
    singleKDePFQ(k,     clauses, vars, percentile,  folder)
    print("KT3DPFQ")
    singleKT3DPFQ(k,    clauses, vars, percentile,  folder)

    for dest in list(glob.glob("GraphEval/*CMCQ*.gml.bz2")):
        G = nx.read_gml(dest, destringizer=int)
        print(MultiGraphToLatexTikz(G, 10, False), file=open(dest.replace("GraphEval", "pics").replace(".gml.bz2", ".tex"), 'w'))

    for dest in glob.glob("GraphEval/*PFQ*.gml.bz2"):
        G = nx.read_gml(dest, destringizer=int)
        print(MultiGraphToLatexTikz(G, 10, clique_sort), file=open(dest.replace("GraphEval", "pics").replace(".gml.bz2", ".tex"), 'w'))


def singleKT3CMCQ(k, c, v, folder="Experiments"):
    """
        k-Sat -Textbook-> 3-SAT -Choi(MIS)-> QUBO
        Input:  k in k-SAT
                v: variable pool
                c: number of clauses
        Output: Metrics as dict
    """
    t1 = time.time()
    sat = createSATinstance(k,c,v)
    t2 = time.time()
    res = {}
    res["Path"] = "KT3CMCQ"
    res["kSAT_k"] = k
    res["kSAT_v"] = v
    res["kSAT_c"] = c
    res["kSAT_actual_v"] = len(getSATVars(sat))
    res["kSAT_actual_lit"] = len(getSATliterals(sat))
    res["kSAT_clauses_per_variable"] = c / len(getSATVars(sat))
    res["kSAT_clauses_per_literal"] = c / len(getSATliterals(sat))

    t3 = time.time()
    sat3 = textbook3SatfromkSat(sat)
    t4 = time.time()
    res["3SAT_c"] = len(sat3)
    res["3SAT_actual_v"] = len(getSATVars(sat3))
    res["3SAT_actual_lit"] = len(getSATliterals(sat3))
    res["3SAT_clauses_per_variable"] = len(sat3) / len(getSATVars(sat3))
    res["3SAT_clauses_per_literal"] = len(sat3) / len(getSATliterals(sat3))

    t5 = time.time()
    Gmis, qubo = ChoisReduction(sat3)
    t6 = time.time()

    print(list(nx.enumerate_all_cliques(Gmis)))
    print(len(sat3))

    pBqubo = PolyBinary(QuboMatToPolyDict(qubo))
    densities = getDensities(pBqubo)
    res["qubo_density_deg1"] = densities[1]
    res["qubo_density_deg2"] = densities[2]
    res["qubo_variables"] = len(qubo)
    res["qubo_monomials"] = len(pBqubo.keys())
    res["qubo_monomials_per_variable"] = len(pBqubo.keys()) / len(qubo)

    res["time_createInstance"] = t2 - t1
    res["time_textbook3SatfromkSat"] = t4 - t3
    res["time_ChoisReduc"] = t6 - t5
    res["time_totalPath"] = (t2 - t1) + (t4 - t3) + (t6 - t5) 
    
    path = folder + "/KT3CMCQ_" + str(k) + "k_" + str(c) + "c_" + str(v) + "v"
    np.save(path + "_qubo.npy", pBqubo)

    nx.write_gml(Gmis, path + "_mis.gml.bz2")

    Gsat = primalGraphFromkSAT(sat)
    nx.write_gml(Gsat, path + "_ksat.gml.bz2")

    G3sat = primalGraphFromkSAT(sat3)
    nx.write_gml(G3sat, path + "_3sat.gml.bz2")

    GQ = primalGraphFromPolynomial(pBqubo)
    nx.write_gml(GQ, path + "_qubo.gml.bz2")

    #df = pd.DataFrame(res, index=[0])
    #df.to_csv(path + ".csv")

    return res

def singleKCMCQ(k, c, v, folder="Experiments"):
    """
        k-Sat -DeMorgan-> PUBO -Fast-> QUBO
        Input:  k in k-SAT
                v: variable pool
                c: number of clauses
        Output: Metrics as dict
    """
    t1 = time.time()
    sat = createSATinstance(k,c,v)
    t2 = time.time()
    res = {}
    res["Path"] = "KCMCQ"
    res["kSAT_k"] = k
    res["kSAT_v"] = v
    res["kSAT_c"] = c
    res["kSAT_actual_v"] = len(getSATVars(sat))
    res["kSAT_actual_lit"] = len(getSATliterals(sat))
    res["kSAT_clauses_per_variable"] = c / len(getSATVars(sat))
    res["kSAT_clauses_per_literal"] = c / len(getSATliterals(sat))

    t3 = time.time()
    Gmis, qubo = ChoisReduction(sat)
    t4 = time.time()

    pBqubo = PolyBinary(QuboMatToPolyDict(qubo))
    densities = getDensities(pBqubo)
    res["qubo_density_deg1"] = densities[1]
    res["qubo_density_deg2"] = densities[2]
    res["qubo_variables"] = len(qubo)
    res["qubo_monomials"] = len(pBqubo.keys())
    res["qubo_monomials_per_variable"] = len(pBqubo.keys()) / len(qubo)

    res["time_createInstance"] = t2 - t1
    res["time_ChoisReduc"] = t4 - t3
    res["time_totalPath"] = (t2 - t1) + (t4 - t3)
    
    path = folder + "/KCMCQ_" + str(k) + "k_" + str(c) + "c_" + str(v) + "v"
    np.save(path + "_qubo.npy", pBqubo)

    nx.write_gml(Gmis, path + "_mis.gml.bz2")

    Gsat = primalGraphFromkSAT(sat)
    nx.write_gml(Gsat, path + "_ksat.gml.bz2")

    GQ = primalGraphFromPolynomial(pBqubo)
    nx.write_gml(GQ, path + "_qubo.gml.bz2")

    #df = pd.DataFrame(res, index=[0])
    #df.to_csv(path + ".csv")

    return res

def singleKDePFQ(k, c, v, p, folder="Experiments"):
    """
        k-Sat -DeMorgan-> PUBO -Fast-> QUBO
        Input:  k in k-SAT
                v: variable pool
                c: number of clauses
                p: percentile for fast quadratisation
        Output: Metrics as dict
    """
    t1 = time.time()
    sat = createSATinstance(k,c,v)
    t2 = time.time()
    res = {}
    res["Path"] = "KDePFQ"
    res["kSAT_k"] = k
    res["kSAT_v"] = v
    res["kSAT_c"] = c
    res["qubo_percentile"] = p
    res["kSAT_actual_v"] = len(getSATVars(sat))
    res["kSAT_actual_lit"] = len(getSATliterals(sat))
    res["kSAT_clauses_per_variable"] = c / len(getSATVars(sat))
    res["kSAT_clauses_per_literal"] = c / len(getSATliterals(sat))

    t3 = time.time()
    pubo = PUBOfromSAT(sat)
    t4 = time.time()

    densitiesPubo = getDensities(pubo)

    for d in range(1, len(densitiesPubo.keys())):
        res["pubo_density_deg" + str(d)] = densitiesPubo[d]

    res["pubo_variables"] = len(pubo.variables)
    res["pubo_monomials"] = len(pubo.keys())
    res["pubo_monomials_per_variable"] = len(pubo.keys()) / len(pubo)

    t5 = time.time()
    qubo, penalty = fastPolyQuadratisation(pubo, max_degree=2, selection_quantile=p)
    qubo_res = PolyBinary(qubo) + PolyBinary(penalty)
    t6 = time.time()

    densitiesQubo = getDensities(qubo_res)

    res["qubo_density_deg1"] = densitiesQubo[1]
    res["qubo_density_deg2"] = densitiesQubo[2]
    res["qubo_variables"] = len(qubo_res.variables)
    res["qubo_monomials"] = len(qubo_res.keys())
    res["qubo_monomials_per_variable"] = len(qubo_res.keys()) / len(qubo_res)


    res["time_createInstance"] = t2 - t1
    res["time_PubofromSat"] = t4 - t3
    res["time_fastQuadratisation"] = t6 - t5
    res["time_totalPath"] = (t2 - t1) + (t4 - t3) + (t6 - t5)
    
    path = folder + "/KDePFQ_" + str(k) + "k_" + str(c) + "c_" + str(v) + "v_" + str(p) + "p"
    np.save(path + "_qubo.npy", qubo_res)

    Gsat = primalGraphFromkSAT(sat)
    nx.write_gml(Gsat, path + "_ksat.gml.bz2")
    
    GP = primalGraphFromPolynomial(pubo, Gsat)
    nx.write_gml(GP, path + "_pubo.gml.bz2")

    GQ = primalGraphFromPolynomial(qubo_res, GP)
    nx.write_gml(GQ, path + "_qubo.gml.bz2")

    #df = pd.DataFrame(res, index=[0])
    #df.to_csv(path + ".csv")

    return res

def singleKT3DPFQ(k, c, v, p, folder="Experiments"):
    """
        k-Sat -Textbook-> 3-SAT -Dobrynin-> PUBO -Fast-> QUBO
        Input:  k in k-SAT
                v: variable pool
                c: number of clauses
                p: percentile for fast quadratisation
        Output: Metrics as dict
    """
    t1 = time.time()
    sat = createSATinstance(k,c,v)
    t2 = time.time()
    res = {}
    res["Path"] = "KT3DPFQ"
    res["kSAT_k"] = k
    res["kSAT_v"] = v
    res["kSAT_c"] = c
    res["qubo_percentile"] = p
    res["kSAT_actual_v"] = len(getSATVars(sat))
    res["kSAT_actual_lit"] = len(getSATliterals(sat))
    res["kSAT_clauses_per_variable"] = c / len(getSATVars(sat))
    res["kSAT_clauses_per_literal"] = c / len(getSATliterals(sat))

    t3 = time.time()
    sat3 = textbook3SatfromkSat(sat)
    t4 = time.time()
    res["3SAT_c"] = len(sat3)
    res["3SAT_actual_v"] = len(getSATVars(sat3))
    res["3SAT_actual_lit"] = len(getSATliterals(sat3))
    res["3SAT_clauses_per_variable"] = len(sat3) / len(getSATVars(sat3))
    res["3SAT_clauses_per_literal"] = len(sat3) / len(getSATliterals(sat3))

    t5 = time.time()
    pubo = PUBOfromSAT(sat3)
    t6 = time.time()

    densitiesPubo = getDensities(pubo)

    for d in range(1, len(densitiesPubo.keys())):
        res["pubo_density_deg" + str(d)] = densitiesPubo[d]

    res["pubo_variables"] = len(pubo.variables)
    res["pubo_monomials"] = len(pubo.keys())
    res["pubo_monomials_per_variable"] = len(pubo.keys()) / len(pubo)

    t7 = time.time()
    qubo, penalty = fastPolyQuadratisation(pubo, max_degree=2, selection_quantile=p)
    qubo_res = PolyBinary(qubo) + PolyBinary(penalty)
    t8 = time.time()

    densitiesQubo = getDensities(qubo_res)

    res["qubo_density_deg1"] = densitiesQubo[1]
    res["qubo_density_deg2"] = densitiesQubo[2]
    res["qubo_variables"] = len(qubo_res.variables)
    res["qubo_monomials"] = len(qubo_res.keys())
    res["qubo_monomials_per_variable"] = len(qubo_res.keys()) / len(qubo_res)

    res["time_createInstance"] = t2 - t1
    res["time_textbook3SatfromkSat"] = t4 - t3
    res["time_PubofromSat"] = t6 - t5
    res["time_fastQuadratisation"] = t8 - t7
    res["time_totalPath"] = (t2 - t1) + (t4 - t3) + (t6 - t5) + (t8 - t7)
    
    path = folder + "/KT3DPFQ_" + str(k) + "k_" + str(c) + "c_" + str(v) + "v_" + str(p) + "p"
    np.save(path + "_qubo.npy", qubo_res)

    Gsat = primalGraphFromkSAT(sat)
    nx.write_gml(Gsat, path + "_ksat.gml.bz2")

    G3sat = primalGraphFromkSAT(sat3, Gsat)
    nx.write_gml(G3sat, path + "_3sat.gml.bz2")

    GP = primalGraphFromPolynomial(pubo, G3sat)
    nx.write_gml(GP, path + "_pubo.gml.bz2")

    GQ = primalGraphFromPolynomial(qubo_res, GP)
    nx.write_gml(GQ, path + "_qubo.gml.bz2")
    
    #df = pd.DataFrame(res, index=[0])
    #df.to_csv(path + ".csv")

    return res