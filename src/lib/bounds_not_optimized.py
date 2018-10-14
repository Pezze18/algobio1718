import numpy as np
import networkx as nx
from lib.core import *
import numpy as np
import bottleneck as bottle
import math
#import tensorflow as tf



###########################################################
####### BOUND_MEANS_ITERATIONS_NOT_OPTIMIZED ###############
###########################################################
def pre_bound_means_iterations_not_optimized(self):
    self.matrix = toMatrix(self, self.G.nodes)

    best_vectors = [ [ [np.ones(len(self.samples)) for i in range(0, self.k)] for cont in range(10000)] for i in range(11)]

    # Distanza 0
    for v in self.G.nodes:
        best_vectors[0][v][0] = np.prod(self.matrix[v])

    sorted_vertices=[(v,best_vectors[0][v][0]) for v in self.G.nodes if best_vectors[0][v][0] !=1]#elimino tutti i nodi inutili che hanno tutte le celle a 1
    sorted_vertices=sorted(sorted_vertices, key= lambda x:x[1])

    sorted_prod=[ x[1] for x in sorted_vertices]
    print(sorted_prod)
    cumsum=np.cumsum(sorted_prod)
    cumsum = cumsum / cumsum[len(cumsum) - 1]

    tresholds = [0,0.000001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    contatori=[0 for i in range(len(tresholds))]
    len_cont = len(contatori)

    j=0
    for i in range(len(cumsum)):
        if cumsum[i]>tresholds[j]:
            contatori[j]=i
            j=j+1
            if(j==len_cont):
                break
    print(contatori)

    M = self.G.copy()
    self.index=0
    self.contatori=contatori
    self.sorted_vertices=[ s[0] for s in sorted_vertices]


    for index in range(len(contatori)):
        if index>0:
            self.G.remove_nodes_from(self.sorted_vertices[self.contatori[index - 1]:self.contatori[index]])
        best_vectors[index]=generate_best_vectors(self, self.G)

    self.best_vectors=best_vectors
    self.G=M


def bound_means_iterations_not_optimized(self, C, vecC):
    dist = self.k - len(C)
    length=len(self.samples)

    bestProd=1
    for g in C:
        bestProd*=self.best_vectors[0][g][0]

    minProd=1
    for g in C:
        if minProd > self.best_vectors[self.index][g][dist] :
            minProd=self.best_vectors[self.index][g][dist]

    bestProd*=minProd
    bestProd=math.pow(bestProd, 1.0/length)
    bestProd*=length

    if bestProd > self.best_score:
        return True
    return False

def update_bound_means_iterations_not_optimized(self):
    for i in range(len(self.contatori)):
        if(self.cont>=self.contatori[i]):
            self.index=i



def generate_best_vectors(self,G):
    #Distanza 0
    best_vectors=[[1 for k in range(self.k)] for v in range(10000)]
    for v in G.nodes:
        best_vectors[v][0]=np.prod(self.matrix[v])

    #Distanza 1
    for v in G.nodes:
        min_prod=1
        for u in G.neighbors(v):
            if min_prod>best_vectors[u][1]:
                min_prod=best_vectors[u][1]
        best_vectors[v][1]=min_prod

    #Distanza da 2 a k-1
    for k in range(2,self.k):
        for v in G.nodes:
            min_prod=1
            for u in G.neighbors(v):
                best_prod_u=best_vectors[u][k-1]*best_vectors[u][0]#best product a dist=1 per se stesso
                if min_prod>best_prod_u:
                    min_prod=best_prod_u
            best_vectors[v][k]=min_prod
    return best_vectors



###################################
####### BOUND_MEANS_NOT_OPTIMEZED ###############
###################################
def pre_bound_means_not_optimized(self):
    self.matrix = toMatrix(self, self.G.nodes)
    self.sorted_vertices = ordinamentoVertici_nobound_migliorato(self)

    #Distanza 0
    best_vectors=[[1 for k in range(self.k)] for v in range(10000)]
    for v in self.G.nodes:
        best_vectors[v][0]=np.prod(self.matrix[v])

    #Distanza 1
    for v in self.G.nodes:
        min_prod=1
        for u in self.G.neighbors(v):
            if min_prod>best_vectors[u][0]:
                min_prod=best_vectors[u][0]
        best_vectors[v][1]=min_prod

    #Distanza da 2 a k-1
    for k in range(2,self.k):
        for v in self.G.nodes:
            min_prod=1
            for u in self.G.neighbors(v):
                best_prod_u=best_vectors[u][k-1]*best_vectors[u][0]#best product a dist=1 per se stesso
                if min_prod>best_prod_u:
                    min_prod=best_prod_u
            best_vectors[v][k]=min_prod

    self.best_vectors=best_vectors

def bound_means_not_optimized(self, C, vecC):
    dist = self.k - len(C)
    length=len(self.samples)

    bestProd=1
    for g in C:
        bestProd*=self.best_vectors[g][0]

    minProd=1
    for g in C:
        if minProd > self.best_vectors[g][dist] :
            minProd=self.best_vectors[g][dist]


    bestProd*=minProd
    bestProd=math.pow(bestProd, 1.0/length)
    bestProd*=length

    if bestProd > self.best_score:
        return True
    return False

def update_bound_means_not_optimized(self):
    return True

######################################################
####### BOUND_KANTOROVICH_NOT_OPTIMIZED ###############
#######################################################
def pre_bound_kantorovich_not_optimized(self):
    self.matrix = toMatrix(self, self.G.nodes)
    for g in self.genes:
        self.matrix[g]=self.matrix[g]
    self.normL2=[len(self.samples) for i in range(10000)]
    for g in self.genes:
        self.normL2[g]=np.linalg.norm(self.matrix[g])
    self.mins=[ [[1 for i in range(self.k)] for j in range(10000)] for index in range(2)]

    M=self.G.copy()
    contatori=[0,156]
    for index in range(2):
        if(index>0):
            self.G.remove_nodes_from(self.sorted_vertices[0:contatori[index]])
        if(index==0):
            num_zeros=0
            self.sorted_vertices=[]
            sorted_vertices_zeros=[]
            sorted_vertices_diff_ones=[]
            for g in self.genes:
                self.mins[index][g][0] = np.min(self.matrix[g])
                if self.mins[index][g][0]==0:
                    num_zeros+=1
                    sorted_vertices_zeros.append(g)
                elif self.mins[index][g][0]!=1:
                    sorted_vertices_diff_ones.append(g)
            print("Numero di vettori che hanno minimo a 0: "+str(num_zeros))
            #Creo Ordinamento
            self.sorted_vertices+=sorted_vertices_zeros
            self.sorted_vertices +=sorted_vertices_diff_ones

            strs=[self.id_to_str[s] for s in sorted_vertices_zeros]
            print(strs)
            print()
            #questo ordinamento mi sarà utile quando terrò conto dell'indice

        num_zeros=0
        for g in self.G.nodes:
            neighbors=self.G.neighbors(g)
            for u in neighbors:
                if self.mins[index][g][1] > self.mins[0][u][0]:
                    self.mins[index][g][1]=self.mins[0][u][0]
            if self.mins[index][g][1]==0:
                num_zeros+=1
        print("Numero di vettori che hanno vicino minimo a 0: " + str(num_zeros))

        #Nota: il max sarà sempre a 1 quindi non è necessario salavarlo
        #$x^T \cdot y \geq min(x) \cdot \min(y) \cdot |x| \cdot |y|$
    self.G=M
    self.index=0
    print("Len sorted vertices: "+str(len(self.sorted_vertices)))

def bound_kantorovich_not_optimized(self,C,vecC):
    dist=self.k-len(C)
    minx=np.min(vecC)
    maxx=np.max(vecC)
    norm_x=np.linalg.norm(vecC)
    if(dist==1):
        miny=1
        min_norm_y=len(self.samples)
        for c in C:
            if miny > self.mins[self.index][c][1]:
                miny=self.mins[self.index][c][1]
            if self.normL2[c] < min_norm_y:
                min_norm_y = self.normL2[c]

        rapy=miny/1
        rapx=minx/maxx
        bestS=rapx*rapy*norm_x*min_norm_y
        #print(bestS)
        if bestS>self.best_score:
            return True
    return False

def update_bound_kantorovich_not_optimized(self):
    if(self.cont>=157):
        self.index=1
    return True

#################################################
####### BOUND_ORDER_NOT_OPTIMIZED ###############
#################################################
def pre_bound_order_not_optimized(self):

    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order(self)
    import pickle
    fileObject=open(self.parameters["bestVectors"],"rb")
    print("bestVectors loading")
    self.best_vectors=pickle.load(fileObject)
    self.max_counts=pickle.load(fileObject)
    fileObject.close()
    print("Fine loading bestVectors")



def bound_order_not_optimized(self,C,vecC):
    dist=self.k-len(C)
    lista=which_diff(vecC)
    dec=np.asarray(lista)
    bests=[]
    if(dist==1 or dist==2):
        for v in C:
            best_vector=self.best_vectors[v][dist]#già ordinato in ordine crescente
            if((len(dec)+len(best_vector) >len(self.samples))):
                #print("Attenzione ! Serve intersezione !")
                dec=np.sort(dec)
                dec=dec[::-1]#ordino vecC(solo valori diversi da 1) in ordine decrescente
                if(len(best_vector)>len(dec)):
                    inters=len(dec)- ( len(self.samples)-len(best_vector) )
                else:
                    inters = len(best_vector) - (len(self.samples) - len(dec))
                #print(len(best_vector))
                bestS=np.sum(best_vector[:len(best_vector)-inters])+np.sum(dec[inters:])
                bestS+=np.dot(best_vector[len(best_vector)-inters:], dec[0:inters])
                #le 3 righe sotto sono equivalenti alle 2 di sopra per il calcolo del bestS ma sono più leggibili
                #dec_reduce=dec[0:len(dec)-inters]
                #best_vector_reduce=best_vector[inters:len(best_vector)]#return a view non una copia
                #bestS=np.sum(dec_reduce)+np.sum(best_vector_reduce)+np.dot(dec[len(dec)-inters:len(dec)], best_vector_reduce[0:inters] )
            else:
                bestS= np.sum(dec)+np.sum(best_vector) + (len(self.samples)-len(dec)-len(best_vector))
            bests.append(bestS)
        bestS=np.min(bests)
        if(bestS > self.best_score):
            return True
    return False

def update_bound_order_not_optimized(self):
    return True
def ordinamentoVertici_bound_order_not_optimized(self):
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
    self.max_counts=[[0 for i in range(self.k)] for i in range(10000)]#quelli che non esistono o che hanno 0 numeri diversi da 1 sono settati a 0,
                                            #ovvero hanno 0 numeri diversi da 1
    for c in cont:
        self.max_counts[c[0]][0]=c[1]

    cont=sorted(cont, key=lambda x: x[1],reverse=True)
    #print(cont)  # guarda distribuzione max_count
    for i in range(len(cont)):
            sorted_vertices.append(cont[i][0])

    #Per ora in standby
    freq = [c[1] for c in cont]
    #Debugging part
    print("Numero di geni con almeno una cella diversa da 1: "+str(len(cont)))
    unique, counts = np.unique(freq, return_counts=True)
    print(np.asarray((unique, counts)).T)

    """
    cumsum=np.cumsum(freq)
    #print(cumsum[len(cumsum) - 1])
    cumsum=cumsum/cumsum[len(cumsum)-1]
    #print(cumsum)
    tresholds=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    contatori = [0 for i in range(len(tresholds))]
    len_cont=len(contatori)

    j=0
    for i in range(len(cumsum)):
        # elimino tutti i nodi da 0 a cumsum[i](quindi cumsum[i] risulta non compreso) e ricaclcolo i bounds
        if cumsum[i]>=tresholds[j]:
            contatori[j]=i
            print()
            j=j+1
            if(j==len_cont):
                break

    #print()
    #print(contatori)
    #print()
    self.contatori=contatori
    self.index=0
    """
    self.sorted_vertices = sorted_vertices
    self.orderedMatrix = orderedMatrix
    print("Fine Ordinamento")


################################################################
####### BOUND_MIN MIGLIORATO_ITERATIONS_NOT_OPTIMIZED ##########
################################################################
def pre_bound_min_migliorato_iterations_not_optimized(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_min(self)#Creo ordinamento(Stesso di bound_min)
    print("Inizio Ordinamento")

    #best_vectors[v][k] dove k indica la dist e 10000 è l'indice massimo per un nodo

    iterations=5
    step=int(len(self.sorted_vertices)/iterations)
    M = self.G.copy()

    best_vectors = [ [ [np.ones(len(self.samples)) for i in range(0, self.k)] for cont in range(10000)] for i in range(iterations)]

    #tresholds=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    tresholds=[i*(1/iterations) for i in range(iterations)]
    self.contatori=[i*step for i in range(iterations)]
    print(self.contatori)

    for i in range(iterations):
        if(i>0):
            self.G.remove_nodes_from(self.sorted_vertices[self.contatori[i-1]:self.contatori[i]])
        for k in range(1, self.k):
            for v in self.G:
                lista_pazienti = [[] for x in range(len(self.samples))]#a ogni sample associo una lista
                for g in self.G.neighbors(v):
                    if k==1:
                        vec = self.matrix[g]
                    else:
                        vec = best_vectors[i][g][k - 1]
                    for j in range(len((lista_pazienti))):
                        if (vec[j] != 1):
                            lista_pazienti[j].append(vec[j])
                for j in range(len(lista_pazienti)):
                    if (len(lista_pazienti[j]) == 0):
                        lista_pazienti[j] = 1
                    else:
                        lista_pazienti[j]=np.asarray(lista_pazienti[j])
                        lista_pazienti[j]=np.min(lista_pazienti[j])# return min_value
                if(k==1):
                    best_vectors[i][v][k] = np.asarray(lista_pazienti)
                else:
                    best_vectors[i][v][k] = np.multiply(best_vectors[i][v][1], np.asarray(lista_pazienti))
        self.best_vectors = best_vectors
    self.G = M
    self.index=0
    print("Fine Ordinamento")

def bound_min_migliorato_iterations_not_optimized(self, C, vecC):  # DA COMPLETARE
    dist = self.k - len(C)
    bests=[]
    for g in C:
        bestS = np.sum(np.multiply(self.best_vectors[self.index][g][dist], vecC))  # self.best_vectors[g][dist]
        bests.append(bestS)
    bestS=np.min(bests)
    if bestS > self.best_score:
         return True
    return False

def update_bound_min_migliorato_iterations_not_optimized(self):
    for i in range(len(self.contatori)):
        if(self.cont>=self.contatori[i]):
            self.index=i

#######################################################
####### BOUND_MIN MIGLIORATO_NOT_OPTIMIZED ############
#######################################################
def pre_bound_min_migliorato_not_optimized(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_min(self)#Creo ordinamento(Stesso di bound_min)
    print("Inizio Ordinamento")
    best_vectors = [[np.ones(len(self.samples)) for i in range(0, self.k)] for cont in range(10000)]
    #best_vectors[v][k] dove k indica la dist e 10000 è l'indice massimo per un nodo

    for v in self.G:
        best_vectors[v] = [x for x in range(self.k)]  # 0...k-1
    for k in range(1, self.k):
        for v in self.G:
            lista_pazienti = [[] for x in range(len(self.samples))]#a ogni sample associo una lista
            for g in self.G.neighbors(v):
                if k==1:
                    vec = self.matrix[g]
                else:
                    vec = best_vectors[g][k - 1]
                for j in range(len((lista_pazienti))):
                    if (vec[j] != 1):
                        lista_pazienti[j].append(vec[j])
            for j in range(len(lista_pazienti)):
                if (len(lista_pazienti[j]) == 0):
                    lista_pazienti[j] = 1
                else:
                    lista_pazienti[j]=np.asarray(lista_pazienti[j])
                    lista_pazienti[j]=np.min(lista_pazienti[j])# return min_value
            if(k==1):
                best_vectors[v][k] = np.asarray(lista_pazienti)
            else:
                best_vectors[v][k] = np.multiply(best_vectors[v][1], np.asarray(lista_pazienti))
    self.best_vectors = best_vectors
    print("Fine Ordinamento")

def bound_min_migliorato_not_optimized(self, C, vecC):  # DA COMPLETARE
    dist = self.k - len(C)
    minBestS=10000
    for g in C:
        bestS = np.sum(
            np.multiply(self.best_vectors[g][dist], vecC))  # self.best_vectors[g][dist]
        if minBestS > bestS:
            minBestS = bestS
        """
        Debugging part
        print(self.best_vectors[g][dist])
        print(vectorization_solution(self,C))
        print("bestS: "+ str(bestS))
        print("ActualScore: "+str(np.dot(vectorization_solution(self,C))))
        print()
        print("bestS: "+str(bestS))
        print("BestScore: "+str(self.best_score))
        """
    if minBestS > self.best_score:
        return True
    return False

def update_bound_min_migliorato_not_optimized(self):
    return True