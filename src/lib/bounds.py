import numpy as np
import networkx as nx
from lib.core import *
import numpy as np
import bottleneck as bottle
import math
from lib.auxiliary_functions import *
#import tensorflow as tf








########################################################
###########BEST VECTORS DISTANZA 1######################
########################################################

def pre_creaBestVectorsDistanza1(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)
    self.best_vectors = [[[], []] for j in range(self.max_node + 1)]
    creaBestVectors1(self)

    import pickle
    f=open("BestVectorsDistanza1","wb")
    pickle.dump(self.best_vectors,f)
    f.close()
    raise ValueError

def creaBestVectors1(self):
    M=self.G.copy()
    for v in self.G.nodes:
        self.best_vectors[v][0]=updateBestVector(self,self.G,v)

def updateBestVector(self,G,v):
    b = []
    neighbors=self.G.neighbors(v)
    lista = [self.max_counts[u] for u in neighbors]
    if len(lista)==0:
        max_count=0
    else:
        max_count = np.max(lista)
    neighbors = self.G.neighbors(v)
    for u in neighbors:
        ord = np.sort(self.orderedMatrix[u])
        idx = np.searchsorted(b, ord)
        b = np.insert(b, idx, ord)
    b = which_diff(b)[0:max_count]#len(self.samples)
    return b

########################################################
###########BEST VECTORS DISTANZA 2######################
########################################################
def pre_creaBestVectorsDistanza2(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    import pickle
    f=open("BestVectorsDistanza1","rb")
    self.best_vectors = pickle.load(f)
    f.close()

    if self.onlyCount:
        self.max_counts = [0 for i in range(10000)]
    else:
        f=open("BestVectorsDistanza2_max_counts","rb")
        self.max_counts=pickle.load(f)
        f.close()


def crea_creaBestVectorsDistanza2(self,C, vec):
    if self.onlyCount==False:
        diff = which_diff(vec)
        neighbors = set()
        for c in C:
            neighbors.update(self.G.neighbors(c))
        for u in neighbors:
            b = self.best_vectors[u][1]
            ord = np.sort(diff)
            ord = ord[0:self.max_counts[u]]
            idx = np.searchsorted(b, ord)
            b = np.insert(b, idx, ord)[0:self.max_counts[u]]
            self.best_vectors[u][1] = b
    else:
        diff = which_diff(vec)
        max = len(diff)
        neighbors = set()
        for c in C:
            neighbors.update(self.G.neighbors(c))
        for u in neighbors:
            if self.max_counts[u] < max:
                self.max_counts[u] = max


def save_creaBestVectorsDistanza2(self):
    if self.onlyCount:
        import pickle
        f = open("BestVectorsDistanza2_max_counts", "wb")
        pickle.dump(self.max_counts, f)
        f.close()
    else:
        import pickle
        f=open("BestVectorsDistanza2","wb")
        pickle.dump(self.best_vectors,f)
        f.close()

def update_creaBestVectorsDistanza2(self):
    return True

########################################################
###########BOUND ORDER IMPROVED#########################
########################################################

def ordinamentoVertici_bound_order_improved(self):
    print("Inizio Ordinamento")
    #ritorna self.sorted_vertices e self.max_counts
    cont = [(i,0) for i in range(10000)]#mi pare sia corretto
    sorted_vertices = []
    orderedMatrix=[[] for i in range(10000)]
    for g in self.genes:
        orderedMatrix[g]=list(which_diff(self.matrix[g]))
        cont.append((g, len(orderedMatrix[g])))  # (IdNodo,max_count)
    cont=[c for c in cont if c[1]!=0]

    #Mi salvo per ogni nodo il proprio max_count(a livello 0 quindi)
    self.max_counts=[0 for i in range(10000)]#quelli che non esistono o che hanno 0 numeri diversi da 1 sono settati a 0,
                                            #ovvero hanno 0 numeri diversi da 1
    for c in cont:
        self.max_counts[c[0]]=c[1]
    self.orderedMatrix=orderedMatrix

    cont=sorted(cont, key=lambda x: x[1],reverse=True)
    for i in range(len(cont)):
            sorted_vertices.append(cont[i][0])

    #Per ora in standby
    freq = [c[1] for c in cont]

    percentiles=np.asarray([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    thresholds=percentiles
    indici,values=calculatePercentiles(self, freq, thresholds)
    print("contatori: "+str(indici))

    self.index=0
    self.contatori = indici
    self.sorted_vertices = sorted_vertices
    self.orderedMatrix = orderedMatrix
    print("Fine Ordinamento")

def pre_bound_order_improved(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    import pickle
    f=open("BestVectorsDistanza2","rb")
    self.best_vectors = pickle.load(f)
    f.close()

def bound_order_improved(self,C,vecC):
    dist=self.k-len(C)
    lista=which_diff(vecC)
    dec=np.asarray(lista)
    bests=[]
    if(dist==1 or dist==2):
        v=C[len(C)-1]
        best_vector=self.best_vectors[v][dist-1]#già ordinato in ordine crescente
        if((len(dec)+len(best_vector) >len(self.samples))):
            dec=np.sort(dec)
            dec=dec[::-1]#ordino vecC(solo valori diversi da 1) in ordine decrescente
            if(len(best_vector)>len(dec)):
                inters=len(dec)- ( len(self.samples)-len(best_vector) )
            else:
                inters = len(best_vector) - (len(self.samples) - len(dec))
            bestS=np.sum(best_vector[:len(best_vector)-inters])+np.sum(dec[inters:])
            bestS+=np.dot(best_vector[len(best_vector)-inters:], dec[0:inters])
        else:
            bestS= np.sum(dec)+np.sum(best_vector) + (len(self.samples)-len(dec)-len(best_vector))
        if(bestS > self.best_score):
            return True
    return False


def update_bound_order_improved(self):
    return True






