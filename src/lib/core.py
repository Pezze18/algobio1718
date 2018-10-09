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
        self.parameters = parameters
        self.sample_size = len(self.samples)
        self.prob = parameters["prob"]
        self.bound= parameters["bound"]
        self.best_score=parameters["best_score"]
        self.best_subgraph = []
        self.levels=[0 for i in range(k+1)]
        self.str_to_id=parameters["str_to_id"]
        self.id_to_str = parameters["id_to_str"]
        self.genes=list(self.G.nodes).copy()
        self.levels = [0 for i in range(k + 1)] # 0...k
        self.pre=getattr(bounds,"pre_"+parameters["method"])
        self.pre(self)
        self.levelsVecUse=True#parameters["levelsVec"]

        if self.prob and self.bound:
            self.bound_function = getattr(bounds, parameters["method"])
            self.update_function = getattr(bounds, "update_" + parameters["method"])

        if self.prob:
            if self.levelsVecUse:
                self.levelsVec = [1 for i in range(k+1)] # 0...k
                self.scoring_function = getattr(bounds, "prob_cover_vec")
            else:
                self.scoring_function = getattr(bounds, "prob_cover")
        else:
            self.scoring_function = getattr(bounds, "score_cover")
            self.best_score=0#Ricordati di cambiare < con > !


    def inspectNode(self, C):
        size=len(C)
        if size>self.k:
            # Prune it if it's too big: this will never actually happen, since we prune whenever len(C)==k
            return True
        elif size<self.k:
            self.levels[size] += 1
            if self.levelsVecUse:
                if (size == 1):
                    self.levelsVec[1] = self.matrix[C[0]]
                else:
                    self.levelsVec[size] = np.multiply(self.levelsVec[size - 1], self.matrix[C[size - 1]])
            # Prune it if the bounding function can't reach the current best, but do it only when using the probability version.
            return self.prob and self.bound and self.bound_function(self,C,self.levelsVec[size])
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
            if(self.levelsVecUse):
                if (size == 1):
                    self.levelsVec[1] = self.matrix[C[0]]
                else:
                    self.levelsVec[size] = np.multiply(self.levelsVec[size - 1],
                                                       self.matrix[C[size - 1]])  # aggiorno vettore per il livello size
                # This means we've reached a leaf.We should evaluate it, as it wasn't previously pruned!
                score = self.scoring_function(self,self.levelsVec[size])
            else:
                score= self.scoring_function(self, C)#scoring_function per prob_cover senza levelsVec
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

            if(self.prob and self.bound):
                self.update_function(self)

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






class Combinatorial:
    def __init__(self, G, samples, k, parameters):
        self.G = G
        self.k = k
        self.samples = samples
        self.sample_size = len(self.samples)
        self.prob = parameters["prob"]
        self.bound= parameters["bound"]
        self.best_score=0
        self.best_subgraph = []
        self.delta=parameters["delta"]

        if(self.prob):
            self.matrix = toMatrix(self, self.G.nodes)

    def combinatorial_algorithm(self):
        if(self.prob):
            self.combinatorial_algorithm_prob()
        else:
            self.combinatorial_algorithm_det()

    def combinatorial_algorithm_det(self):
        G=self.G
        delta=self.delta

        #G = delta_removal(G, delta)

        C = {}
        P_C = set()

        for v in G.nodes():
            C_v = {v}
            P_C_v = set_cover(self, C_v)  # no need to compute this \foreach u

            p_v = {u: set(nx.shortest_path(G, v, u)) for u in G.nodes() if u is not v}

            while len(C_v) < self.k:
                maximum = -1
                l_v_max = set()

                for u in G.nodes() - C_v:  # "-" is an overloaded operator, it means 'difference'
                    l_v = p_v[u]

                    if len(l_v | C_v) <= self.k:  # "|" is an overloaded operator, it means 'union'
                        P_v = set_cover(self, l_v)

                        s = len(P_v - P_C_v)/ len(l_v - C_v)
                        if maximum < s:
                            maximum = s
                            l_v_max = l_v
                C_v = C_v | l_v_max
                P_C_v = set_cover(self, C_v)

            if len(P_C_v) > len(P_C):  # if we've found a better solution, update it and let us know
                C = C_v
                P_C = P_C_v
                print("_________________")
                print()
                print("Best solution updated!")
                print("Current C (ids): ", C)
                self.best_score=len(P_C)
                self.best_subgraph=C
                print("Current P_C (cardinality):", len(P_C))
        return C, P_C

    def combinatorial_algorithm_prob(self):
        G=self.G
        delta = self.delta
        #G = delta_removal(G, delta)

        C = []
        score_C = -1
        for v in G.nodes():
            C_v = [v]
            vecC_v = self.matrix[v]
            score_C_v=np.sum(vecC_v)

            BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=True, creaEtichetta=True,creaDepths=False)

            while len(C_v) < self.k:
                maximum_rapporto = -1
                max_set = []#nodi da aggiungere a C_v ovvero l_v(u)\C_v
                score_max=-1#indica lo score massimo ottenuto da C_v U l_v(u)
                vec_max=[]#Indica il nuovo vettore ottimo ottenuto da C_v U l_v(u)
                          #dove l_v(u) è l_v(u) che massimizza il rapporto

                for k in range(2,self.k+1):
                    for s in self.L[k]:
                        newC, father=findAncestor(self,v,s)

                        vec_complementare=np.divide(self.shortestVec[s],self.matrix[father])
                        vec_union=np.multiply(vecC_v,vec_complementare)
                        score_union=np.sum(vec_union)

                        rapporto = ( score_C_v - score_union)/ len(newC)
                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union
                            vec_max = vec_union
                C_v+=max_set
                score_C_v = score_max
                vecC_v=vec_max

            if score_C_v > score_C:  # if we've found a better solution, update it and let us know
                C = C_v
                score_C = score_C_v
                print("_________________")
                print()
                print("Best solution updated!")
                print("Current C (ids): ", C)
                self.best_score=score_C_v
                self.best_subgraph=C_v
                print("Current P_C (cardinality):", score_C)
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



