import numpy as np
import networkx as nx
from lib.core import *
import numpy as np
import bottleneck as bottle
import math
#import tensorflow as tf

###################################
####### BOUND_ORDER ########
###################################
def ordinamentoVertici_bound_order(self):
    cont = []
    sorted_vertices = []
    for g in self.genes:
        cont.append((g, how_many_diff(self.matrix[g])))  # IdNodo:max_count
    cont=[c for c in cont if c[1]!=0]
    cont=sorted(cont, key=lambda x: x[1],reverse=True)
    #print(cont)  # guarda distribuzione max_count
    for i in range(len(cont)):
        if cont[i][1] != 0:
            sorted_vertices.append(cont[i][0])
    #print("Numero di geni con almeno una cella diversa da 1: "+str(len(cont)))
    # Debegging, guardo freq cumulata
    freq = [c[1] for c in cont]
    unique, counts = np.unique(freq, return_counts=True)
    #print(np.asarray((unique, counts)).T)

    """
    Analisi
    freq = np.asarray(freq)
    freq = freq / np.sum(freq)
    cum = 0
    threshold = [i * 0.1 for i in range(11)]
    j = 0
    for i in range(len(freq)):
        cum += freq[i]
        if cum > threshold[j]:
            print("Threshold :" + str(threshold[j]) + "at node: " + str(i))
            j += 1
    """


###################################
####### BOUND_MIN MIGLIORATO ########
###################################
def pre_bound_min_migliorato(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_min(self)#Creo ordinamento(Stesso di bound_min)
    print("Inizio Ordinamento")
    best_vectors = [[np.ones(len(self.samples)) for i in range(0, self.k)] for cont in range(10000)]
    #best_vectors[v][k] dove k indica la dist e 10000 è l'indice massimo per un nodo

    for v in self.G:
        lista_pazienti = [[] for x in range(len(self.samples))]#a ogni sample associo una lista
        for g in self.G.neighbors(v):
            vec = self.matrix[g]
            for j in range(len((lista_pazienti))):
                if (vec[j] != 1):
                    lista_pazienti[j].append(vec[j])
        for j in range(len((lista_pazienti))):
            lista_pazienti[j].sort()
        for j in range(len((lista_pazienti))):
            if (len(lista_pazienti[j]) == 0):
                lista_pazienti[j] = 1
            else:
                lista_pazienti[j] = lista_pazienti[j][0]
        best_vectors[v] = [x for x in range(self.k )]#0...k-1
        best_vectors[v][1] = np.asarray(lista_pazienti)

    for k in range(2, self.k):
        for v in self.G:
            lista_pazienti = [[] for x in range(len(self.samples))]
            for g in self.G.neighbors(v):
                vec = best_vectors[g][k - 1]
                for j in range(len((lista_pazienti))):
                    if (vec[j] != 1):
                        lista_pazienti[j].append(vec[j])
            for j in range(len((lista_pazienti))):
                lista_pazienti[j].sort()
                if (len(lista_pazienti[j]) == 0):
                    lista_pazienti[j] = 1
                else:
                    lista_pazienti[j] = lista_pazienti[j][0]
            best_vectors[v][k] = np.multiply(best_vectors[v][1], np.asarray(lista_pazienti))
    self.best_vectors = best_vectors
    print("Fine Ordinamento")

    """
    Idem ma con BFS
    self.visit = [False for i in range(0, len(self.G.nodes) + 1)]
    cont=0
    index=0
    for v in self.G:
        cont+=1
        if(cont>100):
            index+=1
            print(index*100)
            cont=0
        best_vectors[v]={}
        #L=BFS(self,v)
        T=nx.bfs_tree(self.G, v)
        L = []
        L.append([])
        L[0].append(v)
        i = 0
        while (len(L[i]) != 0 and i < self.k - 1):  # L[k-1] livello è presente, nonchè l'ultimo
            L.append([])
            for g in L[i]:
                for u in T.successors(g):
                    L[i + 1].append(u)
            i = i + 1


        #print(v)
        lista_pazienti = [[] for x in range(len(self.samples))]
        for k in range(1,self.k):
            for i in range(len(L[k])):#guardiamo nodi al L-iesimo livello
                vec=self.matrix[L[k][i]]
                for j in range(len(self.samples)):
                    value=vec[j]
                    if(value!=1):
                        lista_pazienti[j].append(value)
            best_vector= [1 for x in range(len(self.samples))]
            for j in range(len(self.samples)):
                length2=min(self.k, len(lista_pazienti[j]))
                lista_pazienti[j].sort()
                lista_pazienti[j]=lista_pazienti[j][0:length2]
                length = min(k, len(lista_pazienti[j]))
                min_value=1
                for i in range(length):
                    min_value*=lista_pazienti[j][i]
                best_vector[j]=min_value
            best_vectors[v][k]=np.asarray(best_vector)
    self.best_vectors=best_vectors
    print("Fine Ordinamento")
    """

    """
    #Equivale alla versione base di bound_min ma senza la dinamicità non serve a nulla
    best_vectors=[x for x in range(len(self.samples))]
    for j in range(len(self.samples)):#per ogni riga dove la riga i-esima indice il paziente i-esimo
        lista=[]
        for i in self.matrix:#per ogni gene nel paziente
            if(self.matrix[i][j]!=1):
                lista.append(self.matrix[i][j])
        lista.sort()
        lista=lista[0:self.k]#seleziono i best k
        best_vectors[j]=[1 for x in range(0,self.k)]#0,1,...k-1
        for k in range(1,self.k):
            prod=1
            length=min(k, len(lista) )
            for t in range(length):
                prod*=lista[t]
            best_vectors[j][k]=prod#best value per i-esima componente/paziente in k nodi
            #print(best_vectors[j][k])
    self.best_vectors={}
    for k in range(1, self.k ):
        self.best_vectors[k]=[]
        for j in range(len(self.samples)):
            #print(best_vectors[j])
            #print(self.best_vectors[k])
            self.best_vectors[k].append(best_vectors[j][k])
        self.best_vectors[k]=np.asarray(self.best_vectors[k])
    print( self.best_vectors)
    print("Fine Ordinamento")
    """

def bound_min_migliorato(self, C, vecC):  # DA COMPLETARE
    dist = self.k - len(C)
    bests = []
    for g in C:
        bestS = np.sum(
            np.multiply(self.best_vectors[g][dist], vecC))  # self.best_vectors[g][dist]
        bests.append(bestS)
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
    bestS =np.min(bests)
    if bestS > self.best_score:  # and abs(bestS-self.best_score)>0.00001
        return True
    return False

###################################
########### BOUND_MIN #############
###################################
def pre_bound_min(self):
    # Preparazione delle strutture necessarie
    self.matrix = toMatrix(self, self.G.nodes)
    qs = {
        s: {gene: 1 - self.samples[s][gene]
            for gene in self.samples[s]}
        for s in self.samples
    }
    sorted_qs = {
        s: sorted(qs[s], key=qs[s].get)
        for s in qs
    }

    sort_by_degree = [[t[0], self.how_many_mins(sorted_qs,t[0]), t[1]]
                           for t in sorted(self.G.degree,
                                           key=lambda t: t[1])]

    genes_dict = {
        a[1]: [
            t[0]
            for t in sort_by_degree
            if t[1] == a[1]
        ]
        for a in sort_by_degree
        if a[1] != 0
    }

    """
    genes=set()
    for s in sorted_qs:
        genes=genes.union(sorted_qs[s])
    print("Numero di geni: "+str(len(genes)))
    """
    v_howmany = {t[0]: t[1] for t in sort_by_degree}

    self.sorted_vertices=[]


    #Itero su tutti i nodi
    cont=0
    print("Inizio Ordinamento")
    best_vectors=[[ np.ones(len(self.samples)) for i in range(0,self.k)] for cont in range(len(self.G.nodes))]
    #best_vectors[cont][k] dove k indica la dist
    for cont in range(len(self.G.nodes)):
        #Seleziono nodo
        max_num_min = max(genes_dict)
        i=max_num_min
        #print("i:"+str(i))
        v = genes_dict[i][0]
        self.sorted_vertices.append(v)

        #Creazione dei best_vectors
        for k in range(1,self.k):
            best_vectors[cont][k]=np.ones(len(self.samples))
            cont_j=0
            for j in self.samples:
                prod = 1.0
                length = min(len(sorted_qs[j]),k)
                for i in range(length):
                    prod *= 1.0-self.samples[j][sorted_qs[j][i]]
                best_vectors[cont][k][cont_j]=prod
                cont_j+=1

        #Aggiornamento delle strutture dopo aver eliminato il nodo
        i=max_num_min
        genes_dict[i].remove(v)
        if len(genes_dict[i]) == 0:
            del genes_dict[i]

        for s in sorted_qs:
            try:
                sorted_qs[s].remove(v)
            except ValueError:
                pass

            if len(sorted_qs[s]) > 0:
                new_min = sorted_qs[s][0]
                old_counter = v_howmany[new_min]
                new_counter = old_counter + 1
                try:
                    genes_dict[old_counter].remove(new_min)
                except KeyError:
                    pass

                genes_dict.setdefault(new_counter, [])
                genes_dict[new_counter].append(new_min)
                v_howmany[new_min] = new_counter

                if old_counter > 0 and len(genes_dict[old_counter]) == 0:
                    del genes_dict[old_counter]
        if len(genes_dict)==0:
            break#non è più necessario proseguire
    self.best_vectors=best_vectors
    print("Fine Ordinamento")

def bound_min(self,C):
    dist=self.k-len(C)
    bestS=np.dot(self.best_vectors[self.cont][dist],vectorization_solution(self,C))
    if(bestS>self.best_score):
        return True #prune

def ordinamentoVertici_bound_min(self):
    # Preparazione delle strutture necessarie
    self.matrix = toMatrix(self, self.G.nodes)
    qs = {
        s: {gene: 1 - self.samples[s][gene]
            for gene in self.samples[s]}
        for s in self.samples
    }
    sorted_qs = {
        s: sorted(qs[s], key=qs[s].get)
        for s in qs
    }

    sort_by_degree = [[t[0], self.how_many_mins(sorted_qs, t[0]), t[1]]
                      for t in sorted(self.G.degree,
                                      key=lambda t: t[1])]

    genes_dict = {
        a[1]: [
            t[0]
            for t in sort_by_degree
            if t[1] == a[1]
        ]
        for a in sort_by_degree
        if a[1] != 0
    }

    v_howmany = {t[0]: t[1] for t in sort_by_degree}

    self.sorted_vertices = []
    # Itero su tutti i nodi
    cont = 0
    print("Inizio Ordinamento")
    best_vectors = [[np.ones(len(self.samples)) for i in range(0, self.k)] for cont in range(len(self.G.nodes))]
    for cont in range(len(self.G.nodes)):
        # Seleziono nodo
        max_num_min = max(genes_dict)
        i = max_num_min
        # print("i:"+str(i))
        v = genes_dict[i][0]
        self.sorted_vertices.append(v)

        # Aggiornamento delle strutture dopo aver eliminato il nodo
        i = max_num_min
        genes_dict[i].remove(v)
        if len(genes_dict[i]) == 0:
            del genes_dict[i]

        for s in sorted_qs:
            try:
                sorted_qs[s].remove(v)
            except ValueError:
                pass

            if len(sorted_qs[s]) > 0:
                new_min = sorted_qs[s][0]
                old_counter = v_howmany[new_min]
                new_counter = old_counter + 1
                try:
                    genes_dict[old_counter].remove(new_min)
                except KeyError:
                    pass

                genes_dict.setdefault(new_counter, [])
                genes_dict[new_counter].append(new_min)
                v_howmany[new_min] = new_counter

                if old_counter > 0 and len(genes_dict[old_counter]) == 0:
                    del genes_dict[old_counter]
        if len(genes_dict) == 0:
            break  # non è più necessario proseguire
    print("Fine Ordinamento")

###################################
########### NO_BOUND #############
###################################
def pre_nobound(self):
    self.matrix = toMatrix(self, self.G.nodes)
    #Ordinamento Vertici
    self.sorted_vertices= [t[0] for t in sorted(self.G.degree, key=lambda t: t[1])]
    #self.sorted_vertices= ordinamentoVertici_bound_min(self)  #permette di essere pù veloce perchè salta 40% dei geni
    # self.sorted_vertices= ordinamentoVertici_nobound_migliorato(self)  #permette di essere pù veloce perchè salta 40% dei geni

def ordinamentoVertici_nobound_migliorato(self):
 #unica differenza da ordinamentoVertici_bound_min_migliorato è che non ordino per max_count
    cont=[]
    sorted_vertices=[]
    for g in self.genes:
        cont.append( (g, len(how_many_diff(self.matrix[g]))) ) #IdNodo:max_count
    for i in range(len(cont)):
        if cont[i][1]!=0:
            sorted_vertices.append(cont[i][0])

###################################
########### DET #############
###################################
def pre_det(self):
    #Ordinamento Vertici
    self.sorted_vertices = list(self.G.nodes).copy()

####################################
#########SCORING FUNCTIONS##########
####################################
def toMatrix(self,genes):#matrice in formato qij
    samples=self.samples
    matrix = [[1 for j in range(len(self.samples))] for i in range(10000)]  # 10000 is number maximum of nodes
    for g in genes:
        index=0
        for s in samples:
            if g in samples[s]:
                matrix[g][index]=1-samples[s][g]
            index+=1
        matrix[g]=np.asarray(matrix[g])
    """
    print("Numero di geni: "+str(len(genes)))
    cont=[]
    for g in genes:
        cont2=0
        for i in range(len(self.samples)):
            if(matrix[g][i]==1):
                cont2+=1
        cont.append(cont2)

    contatore=0
    for c in cont:
        if c!=len(self.samples):
            contatore+=1
    print("Numero di geni con valori diversi da 1: "+str(contatore))
    """
    return matrix

def vectorization_solution(self,C):#nunmpy
    v = self.matrix[C[0]]
    for i in range(1, len(C)):
        v = np.multiply(v, self.matrix[C[i]])
    return v

def prob_cover(self,C):#numpy #min version
    v=vectorization_solution(self,C)
    som = np.sum(v)
    return som

def prob_cover_vec(self,vecC):#numpy #min version
    som = np.sum(vecC)
    return som

def prob_cover_old(self,C):#version min
    som=0.0
    samples=self.samples
    for j in samples:
        prod=1.0
        for i in C:
            if i in samples[j]:
                prod*=1.0-samples[j][i]
        som+=prod
    return som

def score_cover(self,C):
    samples = self.samples
    return len([
        p for p in samples
        if len(samples[p].intersection(C))>0
    ])

def set_cover(self,C):
    samples = self.samples
    return set([
        p for p in samples
        if len(samples[p].intersection(C)) > 0
    ])



####################################
######### AUXILIAR FUNCTIONS #######
####################################

def BFS_complete(self):
    L=[ [1 for v in range(self.k)] for i in range(len(self.G))]
    visit = [[False for v in range(len(self.G))] for i in range(len(self.G))]
    pred = [[None for v in range(len(self.G))] for i in range(len(self.G))]
    shortestVec = [[None for v in range(len(self.G))] for i in range(len(self.G))]
    for v in self.G:
        L[v][0]=v
        for g in L[v][0]:
            lista= self.G.neighbors(g)
            for u in lista:
                if visit[v][u]==False:
                    visit[v][u] = True
                    L[v][1].append(u)
                    pred[v][u]=g
                    shortestVec[v][u]=np.multiply(shortestVec[v][u], self.matrix[v])
    for k in range(2,self.k):
        for v in self.G:
            for g in L[v][k-1]:
                lista= self.G.neighbors(g)
                for u in lista:
                    if visit[u][v]==False:
                        visit[u][v]=True
                        L[v][k].append(u)
                        pred[u][v]=g

def BFS(self, s):#metodo ausiliario usabili in diversi bound come major o kantorovich
    #print("Inizio BFS")
    G=self.G
    visit=self.visit
    for i in range(0,len(visit)):
        visit[i]=False
    L=[]
    L.append([])
    L[0].append(s)
    visit[s]=True
    i=0

    while(len(L[i])!=0 and i<self.k-1):#L[k-1] livello è presente, nonchè l'ultimo
        L.append([])
        for v in L[i]:
            for u in G.neighbors(v):
                if visit[u]==False:
                    L[i+1].append(u)
                    visit[u] = True
        i=i+1
    #print(len(L))
    #print(L)
    return L

def how_many_diff(L):
    L=np.asarray(L)
    return np.count_nonzero(L != 1)

def which_diff(L):
    L=np.asarray(L)
    return L[L!=1]

def how_many_diff_old(L):# metodo ausiliario di bound_min(deprecated)
    lista = []
    for n in L:
        if (n != 1):
            lista.append(n)
    return lista

