import networkx as nx
import lib.bounds as bounds
from lib.bounds import *
import numpy as np
import math

class BDDE:
    def __init__(self, G, samples, k, parameters):
        self.G = G
        self.k = k
        self.tree = None
        self.root = None
        self.samples = samples
        self.sample_size = len(self.samples)
        self.prob = parameters["prob"]
        self.bound= parameters["bound"]
        self.best_score=parameters["best_score"]
        self.best_subgraph = []
        self.levels=[0 for i in range(k+1)]
        self.str_to_id=parameters["str_to_id"]
        self.id_to_str = parameters["id_to_str"]
        self.genes=list(self.G.nodes).copy()

        self.pre=getattr(bounds,"pre_"+parameters["method"])
        self.pre(self)

        if self.prob and self.bound:
            self.bound_function = getattr(bounds, parameters["method"])

        if self.prob:
            self.scoring_function = getattr(bounds, "prob_cover")
        else:
            self.scoring_function = getattr(bounds, "score_cover")
            self.best_score=0


    def inspectNode(self, C):
        size=len(C)
        if size>self.k:
            # Prune it if it's too big: this will never actually happen, since we prune whenever len(C)==k
            return True
        elif size<self.k:
            self.levels[size] += 1
            # Prune it if the bounding function can't reach the current best, but do it only when using the probability version.
            return self.prob and self.bound and self.bound_function(self,C)
        else:
            self.levels[size] += 1
            """
            Debegging generico
            print()
            print(C)
            print("ActualScore: " + str(score))
            #Calcolo bestS
            #print("Optimal Score: "+str(bestS))
            print("BestScore: " + str( self.best_score))
            print()
            """
            # This means we've reached a leaf.We should evaluate it, as it wasn't previously pruned!
            score = self.scoring_function(self,C)
            if score<self.best_score:
                self.best_score = score
                self.best_subgraph = C
                print("_________________")
                print()
                print("Best solution updated!")
                print("Current C (ids): ", self.best_subgraph)
                print("Current P_C (cardinality):", self.best_score)
            # No need to go further.
            return True

    def how_many_mins(self,sorted_qs, v):
        count=0
        for s in sorted_qs:
            if sorted_qs[s][0] == v:
                count+=1
        return count

    def printInitialMessage(self):
        if(self.prob):
            if(self.bound):
                print("Probabilistic con bound "+self.bound_function.__name__+" BDDE starts now:")
            else:
                print("Probabilistic (without bound) " + " BDDE starts now:")
        else:
            print("Deterministic (without bound) " + " BDDE starts now:")

    def enumeration_algorithm(self):
        self.printInitialMessage()
        self.cont=0

        for v in self.sorted_vertices:
            #print("v: "+str(v))
            self.root=v
            self.tree=nx.DiGraph()
            self.DEPTH([],v,[])
            self.G.remove_node(v)

            # Removing whatever we don't use anymore
            del self.tree
            del self.root
            self.cont+=1

        print("Numero foglie: ")
        print(self.levels[self.k])
        print("Numero nodi interni visitati in totale: ")
        print(sum(self.levels)-self.levels[self.k])
        print("Numero nodi interni visitati per livello: ")
        print(self.levels)

    def BREADTH(self, S,n,U):
        vn=n.data
        if vn in U:
            return None

        S1=S+[vn]
        if self.inspectNode(S1):
            return None

        n1=n
        for nxx in self.getNodesFromBranch(n):
            n2=self.BREADTH(S1,nxx,U)
            if n2 is not None:
                self.tree.add_edge(n1,n2)
        return n1

    def DEPTH(self, S,v,beta):
        S1=S+[v]
        if self.inspectNode(S1):
            return None

        n=Nodo(v)
        beta1=[]

        xn=self.getxn(S,v)
        xn.sort(reverse=True)
        for i in range(0,len(beta),1):
            n1=self.BREADTH(S1,beta[i],xn)
            if n1 is not None:
                self.tree.add_edge(n,n1)
                beta1.append(n1)
                del n1
        for v in xn:
            n1=self.DEPTH(S1,v,beta1)
            if n1 is not None:
                self.tree.add_edge(n,n1)
                beta1.append(n1)
        return n

    def getNodesFromBranch(self, n):
        W=[]
        neighbors=self.tree[n]
        for v in neighbors:
            if self.tree.has_edge(n,v):
                W.append(v)
        return W

    def getxn(self, S,v):
        L=nx.Graph()
        L.add_nodes_from(self.G.neighbors(v))

        if self.root in L:
            L.remove_node(self.root)
        for n in S:
            L.remove_nodes_from(self.G.neighbors(n))
        return list(L.nodes())


class Nodo(object):
    def __init__(self,data):
        self.data=data

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self)








def combinatorial_algorithm(G, k, patients, delta=0.8):
    G = delta_removal(G, delta)

    C = {}
    P_C = set()

    for v in G.nodes():
        C_v = {v}
        P_C_v = set_cover(patients, C_v)  # no need to compute this \foreach u

        p_v = {u: set(nx.shortest_path(G, v, u)) for u in G.nodes() if u is not v}

        while len(C_v) < k:
            maximum = -1
            l_v_max = set()

            for u in G.nodes() - C_v:  # "-" is an overloaded operator, it means 'difference'
                l_v = p_v[u]

                if len(l_v | C_v) <= k:  # "|" is an overloaded operator, it means 'union'
                    P_v = set_cover(patients, l_v)

                    s = len(P_v - P_C_v)/ len(l_v - C_v)
                    if maximum < s:
                        maximum = s
                        l_v_max = l_v
            C_v = C_v | l_v_max
            P_C_v = set_cover(patients, C_v)

        if len(P_C_v) > len(P_C):  # if we've found a better solution, update it and let us know
            C = C_v
            P_C = P_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            print("Current P_C (cardinality):", len(P_C))

    return C, P_C

def prob_combinatorial_algorithm(G, k, patients, delta=0.8):
    G = delta_removal(G, delta)

    C = {}
    P_C = -1

    for v in G.nodes():
        C_v = {v}
        P_C_v = prob_cover(patients, C_v)  # no need to compute this \foreach u

        p_v = {u: set(nx.shortest_path(G, v, u)) for u in G.nodes() if u is not v}

        while len(C_v) < k:
            maximum = -1
            l_v_max = set()

            for u in G.nodes() - C_v:  # "-" is an overloaded operator, it means 'difference'
                l_v = p_v[u]

                if len(l_v | C_v) <= k:  # "|" is an overloaded operator, it means 'union'
                    P_v = prob_cover(patients, l_v | C_v)

                    s = (P_v - P_C_v)/ len(l_v - C_v)
                    if maximum < s:
                        maximum = s
                        l_v_max = l_v
            C_v = C_v | l_v_max
            P_C_v = prob_cover(patients, C_v)

        if P_C_v > P_C:  # if we've found a better solution, update it and let us know
            C = C_v
            P_C = P_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            print("Current P_C (cardinality):", P_C)

    return C, P_C

def delta_removal(G, delta):
    """
    DEPRECATED. We're not even using this version of the problem.
    :param G: input graph
    :param delta: threshold
    :return: the graph G_I obtained from G by removing any edge with weight < delta
    """

    removes = []

    for (v, u) in G.edges():
        if G[v][u]['weight'] < delta:
            removes.append((v, u))

    G.remove_edges_from(removes)

    return G

