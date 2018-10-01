import networkx as nx
import lib.bounds as bounds
from lib.bounds import *
import numpy as np
import math
import bottleneck as bottle

class BDDE:
    def __init__(self, G, samples, k, parameters):
        self.G = G
        self.M=G.copy()
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
        self.levels = [0 for i in range(k + 1)] # 0...k
        self.pre=getattr(bounds,"pre_"+parameters["method"])
        self.pre(self)

        if self.prob and self.bound:
            self.bound_function = getattr(bounds, parameters["method"])
            self.update_function = getattr(bounds, "update_" + parameters["method"])

        if self.prob:
            self.levelsVec = [1 for i in range(k+1)] # 0...k
            self.scoring_function = getattr(bounds, "prob_cover_vec")
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
            if (size == 1):
                self.levelsVec[1] = self.matrix[C[0]]
            else:
                self.levelsVec[size] = np.multiply(self.levelsVec[size - 1],
                                                   self.matrix[C[size - 1]])  # aggiorno vettore per il livello size
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
            if (size == 1):
                self.levelsVec[1] = self.matrix[C[0]]
            else:
                self.levelsVec[size] = np.multiply(self.levelsVec[size - 1],
                                                   self.matrix[C[size - 1]])  # aggiorno vettore per il livello size
            # This means we've reached a leaf.We should evaluate it, as it wasn't previously pruned!
            score = self.scoring_function(self,self.levelsVec[size])

            diff=list(which_diff(self.levelsVec[size]))
            neighbors=set()
            for c in C:
                neighbors|=set(self.M.neighbors(c))#devo usare grafo originale dato che alcuni vicini potrebbero essere spariti
            for u in neighbors:
                self.lista_current[u]+=diff
                self.max_counts[u][2]=max(self.max_counts[u][2],len(diff))

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

        for v in self.M.nodes:
            if(len(self.lista_current[v])==self.max_counts[v][2]):
                lista=np.asarray(self.lista_current[v])
            else:
                lista = bottle.partition(self.lista_current[v],self.max_counts[v][2])
            remains = np.sort(lista)
            self.best_vectors[v][2]=remains

        import pickle
        fileObject=open("bestVectors2",'wb')
        pickle.dump(self.best_vectors, fileObject)
        pickle.dump(self.max_counts, fileObject)
        fileObject.close()


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

        C = {}
        P_C = -1

        for v in G.nodes():
            C_v = {v}
            P_C_v = prob_cover_max_version(self, list(C_v))  # no need to compute this \foreach u

            p_v = {u: set(nx.shortest_path(G, v, u)) for u in G.nodes() if u is not v}

            while len(C_v) < self.k:
                maximum = -1
                l_v_max = set()

                for u in G.nodes() - C_v:  # "-" is an overloaded operator, it means 'difference'
                    l_v = p_v[u]

                    if len(l_v | C_v) <= self.k:  # "|" is an overloaded operator, it means 'union'
                        P_v = prob_cover_max_version(self, list(l_v | C_v))

                        s = (P_v - P_C_v)/ len(l_v - C_v)
                        if maximum < s:
                            maximum = s
                            l_v_max = l_v
                C_v = C_v | l_v_max
                P_C_v = prob_cover_max_version(self, list(C_v))

            if P_C_v > P_C:  # if we've found a better solution, update it and let us know
                C = C_v
                P_C = P_C_v
                print("_________________")
                print()
                print("Best solution updated!")
                print("Current C (ids): ", C)
                self.best_score=P_C
                self.best_subgraph=C
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



