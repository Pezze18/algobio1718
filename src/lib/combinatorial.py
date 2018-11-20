import networkx as nx
import lib.bounds as bounds
from lib.bounds import *
import numpy as np
import math
import bottleneck as bottle


def combinatorial_algorithm_det_oldVersion(self):
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

def combinatorial_algorithm_det_onlyBFSAndNumpy(self):  # Complementary incluso
    G=self.G
    k=self.k
    patients=self.samples
    self.matrix = toMatrix_det(self, self.genes)

    C = []
    score_C = -1000
    for v in G.nodes():
        C_v = [v]
        score_C_v = score_cover(self,C_v)

        #Creo albero per nodo v
        self.creaRami=False
        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=False, creaEtichetta=True,
                          creaDepths=False, creaRami=self.creaRami)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []  # nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max = -1  # indica lo score massimo ottenuto da C_v U l_v(u)

            for k in range(2, self.k + 1):
                # print(k)
                for s in self.L[k]:
                    if (self.labels[s]):
                        continue

                    if(self.creaRami):
                        ramo = self.branches[s]
                    else:
                        ramo=findRamo(self,v,s)

                    union = set(ramo)
                    newC = list(union.difference(C_v))
                    union = newC+ C_v

                    if len(newC) + len(C_v) <= self.k:
                        score_union = score_cover(self, union)
                        rapporto = (score_union - score_C_v) / len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union

            # print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v += max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c] = True

        if score_C_v > score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score = score_C_v
            self.best_subgraph = C_v
            print("Current score(min_version):", score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_det_BFSAndNumpyAndComplement(self):  # Complementary incluso
    G=self.G
    k=self.k
    patients=self.samples
    self.matrix = toMatrix_det(self, self.genes)

    C = []
    score_C = -1000
    for v in G.nodes():
        C_v = [v]
        score_C_v = score_cover(self,C_v)

        #Creo albero per nodo v
        self.creaRami=True
        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=False, creaEtichetta=True,
                          creaDepths=False, creaRami=self.creaRami)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []  # nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max = -1  # indica lo score massimo ottenuto da C_v U l_v(u)

            for k in range(2, self.k + 1):
                # print(k)
                for s in self.L[k]:
                    if (self.labels[s]):
                        continue

                    if(self.creaRami):
                        ramo = self.branches[s]
                    else:
                        ramo=findRamo(self,v,s)

                    union = set(ramo)
                    newC = list(union.difference(C_v))
                    union = newC+ C_v

                    if len(newC) + len(C_v) <= self.k:
                        score_union = score_cover(self, union)
                        rapporto = (score_union - score_C_v) / len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union

            # print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v += max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c] = True

        if score_C_v > score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score = score_C_v
            self.best_subgraph = C_v
            print("Current score(min_version):", score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_det_BFSAndLevelsVec(self):
    G=self.G
    k=self.k
    patients=self.samples
    self.matrix = toMatrix_det(self, self.genes)

    C = []
    score_C = -1000
    for v in G.nodes():#[8821]:#
        C_v = [v]
        vecC_v = self.matrix[v]
        score_C_v= score_cover_vec(self, vecC_v)

        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=True, creaEtichetta=True,creaDepths=False,creaRami=False)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []#nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max=-1#indica lo score massimo ottenuto da C_v U l_v(u)
            vec_max=[]#Indica il nuovo vettore ottimo ottenuto da C_v U l_v(u)
                      #dove l_v(u) è l_v(u) che massimizza il rapporto

            for k in range(2,self.k+1):
                #print(k)
                for s in self.L[k]:
                    if(self.labels[s]):
                        continue
                    #print(s)
                    newC, father=findAncestor(self,v,s)
                    if len(newC)+len(C_v)<=self.k:
                        if father != v:
                            vec_complementare = self.shortestVec[s] - self.shortestVec[father]
                        else:
                            #in BFS_complete_nodo la radice ha vettore associato a tutti 0
                            vec_complementare = self.shortestVec[s]

                        vec_union=vecC_v + vec_complementare
                        score_union=score_cover_vec(self, vec_union)

                        rapporto = ( score_union-score_C_v )/ len(newC)

                        """if(newC[0]==6780):
                            print(self.shortestVec[s])
                            print(self.matrix[6780])
                            print(self.matrix[father])
                            print(vec_complementare)
                            print(score_union)
                            print("____________________-")"""

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union
                            vec_max = vec_union

            #print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v+=max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c]=True
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
            print("Current score(min_version):",  score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_prob_oldVersion(self):
    #G = delta_removal(G, delta)
    G=self.G
    k=self.k
    patients=self.samples

    C = {}
    P_C = -1

    for v in G.nodes():
        C_v = {v}
        P_C_v = prob_cover_old(self, C_v)  # no need to compute this \foreach u

        p_v = {u: set(nx.shortest_path(G, v, u)) for u in G.nodes() if u is not v}

        while len(C_v) < k:
            maximum = -1
            l_v_max = set()

            for u in G.nodes() - C_v:  # "-" is an overloaded operator, it means 'difference'
                l_v = p_v[u]

                if len(l_v | C_v) <= k:  # "|" is an overloaded operator, it means 'union'
                    P_v = prob_cover_old(self, l_v | C_v)

                    s = (P_v - P_C_v)/ len(l_v - C_v)
                    if maximum < s:
                        maximum = s
                        l_v_max = l_v
            C_v = C_v | l_v_max
            P_C_v = prob_cover_old(self, C_v)

        if P_C_v > P_C:  # if we've found a better solution, update it and let us know
            C = C_v
            P_C = P_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            print("Current P_C (cardinality):", P_C)

    return C, P_C

def combinatorial_algorithm_prob_onlyBFS(self):  # Complementary incluso
    G = self.G
    delta = self.delta
    # G = delta_removal(G, delta)

    C = []
    score_C = 1000
    for v in G.nodes():  # [8821]:#
        C_v = [v]
        score_C_v = prob_cover_old(self,C_v)

        #Creo albero per nodo v
        self.creaRami=False
        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=False, creaEtichetta=True,
                          creaDepths=False, creaRami=self.creaRami)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []  # nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max = -1  # indica lo score massimo ottenuto da C_v U l_v(u)

            for k in range(2, self.k + 1):
                # print(k)
                for s in self.L[k]:
                    if (self.labels[s]):
                        continue

                    if(self.creaRami):
                        ramo = self.branches[s]
                    else:
                        ramo=findRamo(self,v,s)

                    union = set(ramo)
                    newC = list(union.difference(C_v))
                    union = newC+ C_v

                    if len(newC) + len(C_v) <= self.k:
                        score_union = prob_cover_old(self, union)
                        rapporto = (score_C_v - score_union) / len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union

            # print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v += max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c] = True

        if score_C_v < score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score = score_C_v
            self.best_subgraph = C_v
            print("Current score(min_version):", score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_prob_onlyBFSAndProbCover(self):#Complementary incluso
    G=self.G
    delta = self.delta
    #G = delta_removal(G, delta)

    C = []
    score_C = 1000
    for v in G.nodes():#[8821]:#
        C_v = [v]
        vecC_v = self.matrix[v]
        score_C_v=np.sum(vecC_v)

        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=True, creaEtichetta=True,creaDepths=False,creaRami=True)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []#nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max=-1#indica lo score massimo ottenuto da C_v U l_v(u)
            vec_max=[]#Indica il nuovo vettore ottimo ottenuto da C_v U l_v(u)
                      #dove l_v(u) è l_v(u) che massimizza il rapporto

            for k in range(2,self.k+1):
                #print(k)
                for s in self.L[k]:
                    if(self.labels[s]):
                        continue

                    union = set(self.branches[s])
                    newC = list(union.difference(C_v))
                    union = newC+ C_v

                    if len(newC)+len(C_v)<=self.k:
                        vec_complementare=vectorization_solution(self,newC)
                        vec_union=np.multiply(vecC_v,vec_complementare)
                        score_union=np.sum(vec_union)

                        rapporto = ( score_C_v - score_union)/ len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union
                            vec_max = vec_union

            #print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v+=max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c]=True
            vecC_v=vec_max

        if score_C_v < score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score=score_C_v
            self.best_subgraph=C_v
            print("Current score(min_version):",  score_C)
    return self.best_subgraph, self.best_score

#sfruttando il fatto che trovo velocemente newC dovrebbe essere più veloce
def combinatorial_algorithm_prob_onlyBFSAndComplement(self):  # Complementary incluso
    G = self.G
    delta = self.delta
    # G = delta_removal(G, delta)

    C = []
    score_C = 1000
    for v in G.nodes():  # [8821]:#
        C_v = [v]
        score_C_v = prob_cover_old(self,C_v)

        #Creo albero per nodo v
        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=False, creaEtichetta=True,
                          creaDepths=False, creaRami=False)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []  # nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max = -1  # indica lo score massimo ottenuto da C_v U l_v(u)

            for k in range(2, self.k + 1):
                # print(k)
                for s in self.L[k]:
                    if (self.labels[s]):
                        continue

                    newC, father = findAncestor(self, v, s)

                    union = newC+ C_v

                    if len(newC) + len(C_v) <= self.k:
                        score_union = prob_cover_old(self, union)
                        rapporto = (score_C_v - score_union) / len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union

            # print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v += max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c] = True

        if score_C_v < score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score = score_C_v
            self.best_subgraph = C_v
            print("Current score(min_version):", score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_prob_BFSAndProbCoverAndComplement(self):#Complementary incluso
    G=self.G
    delta = self.delta
    #G = delta_removal(G, delta)

    C = []
    score_C = 1000
    for v in G.nodes():#[8821]:#
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
                #print(k)
                for s in self.L[k]:
                    if(self.labels[s]):
                        continue
                    #print(s)
                    newC, father=findAncestor(self,v,s)
                    if len(newC)+len(C_v)<=self.k:
                        """
                        if father!=v:
                            vec_complementare=np.divide(self.shortestVec[s],self.matrix[father])
                            vec_complementare[vec_complementare == np.inf]=0
                        else:
                            vec_complementare=self.shortestVec[s]

                        vec_union=np.multiply(vecC_v,vec_complementare)
                        score_union=np.sum(vec_union)
                        """
                        vec_complementare=vectorization_solution(self,newC)
                        vec_union=np.multiply(vecC_v,vec_complementare)
                        score_union=np.sum(vec_union)

                        rapporto = ( score_C_v - score_union)/ len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union
                            vec_max = vec_union

            #print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v+=max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c]=True
            vecC_v=vec_max

        if score_C_v < score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score=score_C_v
            self.best_subgraph=C_v
            print("Current score(min_version):",  score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_prob_BFSAndLevelsVec(self):#Complementary incluso
    G=self.G
    delta = self.delta
    #G = delta_removal(G, delta)

    import warnings
    warnings.filterwarnings("error")

    C = []
    score_C = 1000
    for v in G.nodes():#[8821]:#
        C_v = [v]
        vecC_v = self.matrix[v]
        score_C_v=np.sum(vecC_v)

        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=True, creaEtichetta=True,creaDepths=False,creaRami=False)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []#nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max=-1#indica lo score massimo ottenuto da C_v U l_v(u)
            vec_max=[]#Indica il nuovo vettore ottimo ottenuto da C_v U l_v(u)
                      #dove l_v(u) è l_v(u) che massimizza il rapporto

            for k in range(2,self.k+1):
                #print(k)
                for s in self.L[k]:
                    if(self.labels[s]):
                        continue
                    #print(s)
                    newC, father=findAncestor(self,v,s)
                    if len(newC)+len(C_v)<=self.k:
                        if father != v:
                            try:
                                vec_complementare = np.divide(self.shortestVec[s], self.shortestVec[father])
                            except:
                                zeros_indices = self.shortestVec[father] == 0
                                vec_complementare[zeros_indices]=0
                        else:
                            #in BFS_complete_nodo la radice ha vettore associato a tutti 1
                            vec_complementare = self.shortestVec[s]

                        vec_union=np.multiply(vecC_v,vec_complementare)
                        score_union=np.sum(vec_union)

                        rapporto = ( score_C_v - score_union)/ len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union
                            vec_max = vec_union

            #print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v+=max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c]=True
            vecC_v=vec_max

        if score_C_v < score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score=score_C_v
            self.best_subgraph=C_v
            print("Current score(min_version):",  score_C)
    return self.best_subgraph, self.best_score

def combinatorial_algorithm_prob_BFSAndLevelsVecAndPruning(self):#Complementary incluso
    G=self.G
    delta = self.delta
    #G = delta_removal(G, delta)

    C = []
    score_C = 1000
    for v in G.nodes():#[8821]:#
        C_v = [v]
        vecC_v = self.matrix[v]
        score_C_v=np.sum(vecC_v)

        if self.bound_function(self,C,vecC_v):
            #print("break")
            continue

        BFS_complete_node(self, v, creaLevels=True, creaPredecessori=True, creaVec=True, creaEtichetta=True,creaDepths=False,creaRami=False)

        while len(C_v) < self.k:
            maximum_rapporto = -1
            max_set = []#nodi da aggiungere a C_v ovvero l_v(u)\C_v
            score_max=-1#indica lo score massimo ottenuto da C_v U l_v(u)
            vec_max=[]#Indica il nuovo vettore ottimo ottenuto da C_v U l_v(u)
                      #dove l_v(u) è l_v(u) che massimizza il rapporto

            for k in range(2,self.k+1):
                #print(k)
                for s in self.L[k]:
                    if(self.labels[s]):
                        continue
                    #print(s)
                    newC, father=findAncestor(self,v,s)
                    if len(newC)+len(C_v)<=self.k:
                        zeros_indices=self.shortestVec[s]==0
                        if father != v:
                            vec_complementare = np.divide(self.shortestVec[s], self.matrix[father])
                            vec_complementare[zeros_indices]=0
                        else:
                            #in BFS_complete_nodo la radice ha vettore associato a tutti 1
                            vec_complementare = self.shortestVec[s]

                        vec_union=np.multiply(vecC_v,vec_complementare)
                        score_union=np.sum(vec_union)

                        rapporto = ( score_C_v - score_union)/ len(newC)

                        if maximum_rapporto < rapporto:
                            maximum_rapporto = rapporto
                            max_set = newC
                            score_max = score_union
                            vec_max = vec_union

            #print("Aggiorno C_v: "+str(C_v)+" con: "+str(max_set))
            C_v+=max_set
            score_C_v = score_max
            for c in max_set:
                self.labels[c]=True
            vecC_v=vec_max

            if self.bound_function(self, C_v, vecC_v):
                #print("break")
                break

        if score_C_v < score_C:  # if we've found a better solution, update it and let us know
            C = C_v
            score_C = score_C_v
            print("_________________")
            print()
            print("Best solution updated!")
            print("Current C (ids): ", C)
            self.best_score=score_C_v
            self.best_subgraph=C_v
            print("Current score(min_version):",  score_C)
    return self.best_subgraph, self.best_score

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