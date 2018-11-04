import numpy as np
import networkx as nx
from lib.core import *
import numpy as np
import bottleneck as bottle
import math
from lib.auxiliary_functions import *
#import tensorflow as tf



#################################################################################
###########BEST VECTORS DISTANZE SUPERIORI#######################################
#################################################################################
def pre_creaBestVectorsDistanzeSuperiori(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    #Creo tabella iniziale, ogni paziente ha una lista di tuple,
    #  dove ogni tupla è formata da (gene,valore)
    tabella=[ [] for j in range(len(self.samples))]
    for g in self.genes:
        for j in range(len(self.samples)):
            if  self.matrix[g][j]!=1:
                tabella[j].append((g, self.matrix[g][j]))
    for j in range(len(self.samples)):#gli ordino in modo da prenderli ed eliminarli facilmente
        tabella[j].sort(key=lambda x:x[1])
        #print(tabella[j])

    #Importo best_vectors per distanza 1
    import pickle
    filename=self.parameters["bestVectors"]+"/"+"BestVectorsDistanza"+str(self.k-1)
    f=open(filename,"rb")
    self.best_vectors = pickle.load(f)
    f.close()
    dist=self.k

    self.cont=0
    for v in self.sorted_vertices:
        deleteFromTable(self, tabella, v)
        b=updateBestVectorTable(self, tabella, dist)
        self.best_vectors[0][self.cont][dist-1]=b
        self.cont += 1

    import pickle
    filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k )
    f = open(filename, "wb")
    pickle.dump(self.best_vectors,f)
    f.close()
    raise ValueError





#################################################################################
###########BEST VECTORS DISTANZA 1 ITERATIONS PERCENTILES########################
#################################################################################

def pre_creaBestVectorsDistanza1_iterations_percentiles(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)
    self.best_vectors = [ [ [[] for n in range(self.k)] for j in range(self.max_node + 1)]  for i in range(len(self.contatori))]
    creaBestVectors1_iterations_percentiles(self)

    #print(self.best_vectors[0][4322])
    #print(self.best_vectors[3][4322])

    import pickle
    filename=self.parameters["bestVectors"]+"/"+"BestVectorsDistanza1"
    f=open(filename,"wb")
    pickle.dump(self.best_vectors,f)
    f.close()
    raise ValueError

def creaBestVectors1_iterations_percentiles(self):
    M=self.G.copy()

    #Calcolo thresholds da usare
    percentiles = [i * 0.1 for i in range(0, 10)]
    totale = []
    for g in self.G.nodes:
        totale += list(self.matrix[g][self.matrix[g] < 1])
    totale=np.sort(totale)
    indici, values = calculatePercentiles(self, totale, percentiles)
    values=[0]+values+[1]
    self.thresholds=values
    thresholds=self.thresholds
    print("thresholds: "+str(thresholds))

    #Calcolo di max_counts_percentiles per ogni nodo/vettore
    self.max_counts_percentiles = [ [0 for j in range(len(thresholds)-1)] for i in range(self.max_node+1)]
    for g in self.G.nodes:
        counts = np.histogram(self.orderedMatrix[g], thresholds)[0]
        for i in range(len(thresholds)-1):
            if self.max_counts_percentiles[g][i] < counts[i]:
                self.max_counts_percentiles[g][i] = counts[i]

    for i in range(len(self.orderedMatrix)):
        self.orderedMatrix[i]=np.asarray(self.orderedMatrix[i])

    #Calcolo di best vectors a distanza 1 per index 0
    M=self.G.copy()
    self.index = 0
    for v in self.G:
        self.best_vectors[self.index][v][0] = updateBestVector_iterations_percentiles(self,self.G, v)

    # Calcolo di best vectors a distanza 1 per index > 0
    for index in range(1,len(self.contatori)):
        self.index=index
        neighbors=set()
        for v in self.sorted_vertices[self.contatori[index-1]:self.contatori[index]] :
            neighbors.update(self.G.neighbors(v))

        self.G.remove_nodes_from(self.sorted_vertices[self.contatori[index-1]:self.contatori[index]])
        for v in self.G:
            if v in neighbors:
                self.best_vectors[index][v][0] = updateBestVector_iterations_percentiles(self,self.G, v)
            else:
                self.best_vectors[index][v][0] = self.best_vectors[index-1][v][0]
    self.G=M


def updateBestVector_iterations_percentiles(self,G, v):
    max=0
    max_percentiles=[0 for i in range(len(self.thresholds)-1)]
    ordered_percentiles = [np.asarray([]) for i in range(len(self.thresholds)-1)]

    #Calcolo max_percentiles(per ogni percentile calcolo il max_count) e prendo tutti i valori e gli metto in ordered_percentiles
    for u in G.neighbors(v):
        if max < self.max_counts[u]:
            max=self.max_counts[u]
        indices = np.digitize(self.orderedMatrix[u], self.thresholds) - 1
        for i in range(len(self.thresholds)-1):
            if max_percentiles[i] < self.max_counts_percentiles[u][i]:
                max_percentiles[i] = self.max_counts_percentiles[u][i]
            ordered_percentiles[i]=np.concatenate( [ordered_percentiles[i], self.orderedMatrix[u][indices==i] ])

    #Seleziono i valori più piccoli per ogni ordered_percentiles
    for i in range(len(self.thresholds)-1):
        length=max_percentiles[i]
        if length==0:
            ordered_percentiles[i]=[]
        else:
            ordered_percentiles[i] = bottle.partition(ordered_percentiles[i], length - 1)
            ordered_percentiles[i] = ordered_percentiles[i][0:length]
            ordered_percentiles[i] = np.sort(ordered_percentiles[i])

    #if v==4322:
    #    print("max_counts: "+str(max)+"  index: "+str(self.index))

    remains = max
    index_perc=0
    best_vector=np.asarray([])
    while (remains > 0):
        if remains > len(ordered_percentiles[index_perc]):
            best_vector = np.concatenate((best_vector, ordered_percentiles[index_perc]))
            remains -= len(ordered_percentiles[index_perc])
        else:
            best_vector = np.concatenate((best_vector, ordered_percentiles[index_perc][0:remains]))
            remains=0
        index_perc+=1
    return best_vector


##########################################################################################
###########BEST VECTORS DISTANZA 2 ITERATIONS PERCENTILES SINGOLO#########################
##########################################################################################
def pre_creaBestVectorsDistanza_iterations_percentiles_singolo(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    #Importo best_vectors per distanza 1
    import pickle
    filename=self.parameters["bestVectors"]+"/"+"BestVectorsDistanza"+str(self.k-1)
    f=open(filename,"rb")
    self.best_vectors = pickle.load(f)
    f.close()

    #Calcolo percentili per distanza 2 usando il miglior vettore trovato a distanza 2(euristica)
    best_solutions=[ [], [6780], [6780, 8821], [6780, 8821, 7031], [2622, 6793, 6780, 8821] ]
    best_set=best_solutions[self.k]
    a=np.sort(which_diff(vectorization_solution(self,best_set)))
    percentiles = [i * 0.1 for i in range(0, 10)]
    indici, values = calculatePercentiles(self, a, percentiles)
    values=[0]+values+[1]
    self.thresholds=values
    print("thresholds: " + str(self.thresholds))

    if self.onlyCount:
        self.max_counts =[ 0  for j in range(len(self.contatori))]
        self.max_counts_percentiles = [ [0 for s in range(len(self.thresholds)-1)]  for j in range(len(self.contatori))]
    else:
        #Importo max_counts e max_counts_percentiles per creare best_vectors a distanza 2
        filename=self.parameters["bestVectors"]+"/"+"BestVectorsDistanza"+str(self.k)+"_max_counts"
        f=open(filename,"rb")
        self.max_counts=pickle.load(f)
        self.max_counts_percentiles=pickle.load(f)
        f.close()

        self.ordered_percentiles = [ [ [] for s in range(len(self.thresholds) - 1)] for j in range(len(self.contatori))]
    self.index=0


def crea_creaBestVectorsDistanza_iterations_percentiles_singolo(self,C, vec):
    if self.onlyCount==False:
        diff = which_diff(vec)
        indices = np.digitize(diff, self.thresholds) - 1

        for i in range(len(self.thresholds) - 1):
            b=self.ordered_percentiles[self.index][i]
            if self.max_counts_percentiles[self.index][i]>0:
                ord=np.sort(diff[indices==i])
                ord = ord[0:self.max_counts_percentiles[self.index][i]]
                idx=np.searchsorted(b, ord )
                #Fare il subset solo ogni tanto potrebbe far risparmiare tempo ?
                b=np.insert(b, idx, ord)[0:self.max_counts_percentiles[self.index][i]]
                self.ordered_percentiles[self.index][i] = b
    else:
        diff = which_diff(vec)
        max = len(diff)
        indices = np.digitize(diff, self.thresholds) - 1
        max_percentiles = [0 for i in range(len(self.thresholds) - 1)]
        unique, counts = np.unique(indices, return_counts=True)
        for i in range(len(unique)):
            u = unique[i]
            max_percentiles[u] = counts[i]

        if self.max_counts[self.index] < max:
            self.max_counts[self.index] = max
        for i in range(len(self.thresholds) - 1):
            if max_percentiles[i] > self.max_counts_percentiles[self.index][i]:
                self.max_counts_percentiles[self.index][i] = max_percentiles[i]



def save_creaBestVectorsDistanza_iterations_percentiles_singolo(self):
    if self.onlyCount:
        cont = 0
        for index in range(len(self.contatori) - 2, -1, -1):
            if self.max_counts[index] < self.max_counts[index + 1]:
                self.max_counts[index] = self.max_counts[index + 1]
                cont=cont+1
            for i in range(len(self.thresholds) - 1):
                if  self.max_counts_percentiles[index][i] < self.max_counts_percentiles[index+1][i]:
                    self.max_counts_percentiles[index][i] = self.max_counts_percentiles[index+1][i]
                    cont = cont + 1
        print("cont: "+str(cont))

        azzerati=0
        for index in range(len(self.contatori) - 1, -1, -1):
            remains = self.max_counts[index]
            for i in range(len(self.thresholds) - 1):
                if remains > 0:
                    if remains > self.max_counts_percentiles[index][i]:
                        remains -= self.max_counts_percentiles[index][i]
                    else:
                        self.max_counts_percentiles[index][i] = remains
                        remains = 0
                else:
                    self.max_counts_percentiles[index][i] = 0
                    azzerati += 1
        print("azzerati: "+str(azzerati))

        import pickle
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k) + "_max_counts"
        f = open(filename, "wb")
        pickle.dump(self.max_counts, f)
        pickle.dump(self.max_counts_percentiles, f)
        f.close()
    else:
        #Devo aggiornare ordered_percentiles
        cont=0
        for index in range(len(self.contatori) - 2, -1, -1):
            for i in range(len(self.thresholds) - 1):
                a = self.ordered_percentiles[index][i]
                my_values = self.ordered_percentiles[index+1][i]
                if len(my_values)>0:#continua, altrimenti non c'è niente da aggiungere
                    idx = np.searchsorted(a, my_values)
                    num = np.count_nonzero(idx >= len(a))
                    if num<len(my_values):#allora almeno un valore va inserito
                        cont=cont+1
                        b = np.insert(a, idx, my_values)
                        b = b[0:self.max_counts_percentiles[index][i]]
                        #if (len(b) > 0):#dato che num<len(my_values) quest'altra condizione non servirebbe ma giusto per essere sicuri
                        self.ordered_percentiles[index][i] = b
        print("cont: " + str(cont))

        #Dai ordered_percentiles creo best_vectors
        for index in range(len(self.contatori)):
            b_v= np.asarray([])
            for i in range(len(self.thresholds) - 1):
                if self.max_counts_percentiles[index][i] > 0:
                    b_v = np.concatenate((b_v, self.ordered_percentiles[index][i][
                                               0:self.max_counts_percentiles[index][i]]))
            self.best_vectors[index][0][self.k-1] = b_v
        del self.ordered_percentiles


        import pickle
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k)
        f=open(filename,"wb")
        pickle.dump(self.best_vectors,f)
        f.close()

def update_creaBestVectorsDistanza_iterations_percentiles_singolo(self):
    index=self.index
    for i in range(len(self.contatori)):
        if self.cont>=self.contatori[i]:
            index=i
    if index!=self.index:
        self.index = index
        print(self.index)
    print(self.cont)
    print(self.levels)

###############################################################################
###########BEST VECTORS DISTANZA 2 ITERATIONS PERCENTILES######################
###############################################################################
def pre_creaBestVectorsDistanza_iterations_percentiles(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    #Importo best_vectors per distanza 1
    import pickle
    filename=self.parameters["bestVectors"]+"/"+"BestVectorsDistanza"+str(self.k-1)
    f=open(filename,"rb")
    self.best_vectors = pickle.load(f)
    f.close()

    #Calcolo percentili per distanza 2 usando il miglior vettore trovato a distanza 2(euristica)
    best_set= [6780, 8821]
    a=np.sort(which_diff(vectorization_solution(self,best_set)))
    percentiles = [i * 0.1 for i in range(0, 10)]
    indici, values = calculatePercentiles(self, a, percentiles)
    values=[0]+values+[1]
    self.thresholds=values
    print("thresholds: " + str(self.thresholds))

    if self.onlyCount:
        self.max_counts =[ [0 for i in range(self.max_node+1)]  for j in range(len(self.contatori))]
        self.max_counts_percentiles = [[ [0 for s in range(len(self.thresholds)-1)] for i in range(self.max_node + 1)] for j in range(len(self.contatori))]
    else:
        #Importo max_counts e max_counts_percentiles per creare best_vectors a distanza 2
        filename=self.parameters["bestVectors"]+"/"+"BestVectorsDistanza"+str(self.k)+"_max_counts"
        f=open(filename,"rb")
        self.max_counts=pickle.load(f)
        self.max_counts_percentiles=pickle.load(f)
        f.close()

        self.ordered_percentiles = [[[[] for s in range(len(self.thresholds) - 1)] for i in range(self.max_node + 1)]
                                       for j in range(len(self.contatori))]

    self.index=0


def crea_creaBestVectorsDistanza_iterations_percentiles(self,C, vec):
    if self.onlyCount==False:
        diff = which_diff(vec)
        indices = np.digitize(diff, self.thresholds) - 1

        neighbors = set()
        for c in C:
            neighbors.update(self.G.neighbors(c))

        for u in neighbors:
            for i in range(len(self.thresholds) - 1):
                b=self.ordered_percentiles[self.index][u][i]
                if self.max_counts_percentiles[self.index][u][i]>0:
                    ord=np.sort(diff[indices==i])
                    ord = ord[0:self.max_counts_percentiles[self.index][u][i]]
                    idx=np.searchsorted(b, ord )
                    b=np.insert(b, idx, ord)[0:self.max_counts_percentiles[self.index][u][i]]
                    self.ordered_percentiles[self.index][u][i] = b
    else:
        diff = which_diff(vec)
        max = len(diff)
        indices = np.digitize(diff, self.thresholds) - 1
        max_percentiles = [0 for i in range(len(self.thresholds) - 1)]
        unique, counts = np.unique(indices, return_counts=True)
        for i in range(len(unique)):
            u = unique[i]
            max_percentiles[u] = counts[i]

        neighbors = set()
        for c in C:
            neighbors.update(self.G.neighbors(c))

        for u in neighbors:
            if self.max_counts[self.index][u] < max:
                self.max_counts[self.index][u] = max
            for i in range(len(self.thresholds) - 1):
                if max_percentiles[i] > self.max_counts_percentiles[self.index][u][i]:
                    self.max_counts_percentiles[self.index][u][i] = max_percentiles[i]



def save_creaBestVectorsDistanza_iterations_percentiles(self):
    if self.onlyCount:
        cont = 0
        for index in range(len(self.contatori) - 2, -1, -1):
            for v in self.genes:
                if self.max_counts[index][v] < self.max_counts[index + 1][v]:
                    self.max_counts[index][v] = self.max_counts[index + 1][v]
                    cont=cont+1
                for i in range(len(self.thresholds) - 1):
                    if  self.max_counts_percentiles[index][v][i] < self.max_counts_percentiles[index+1][v][i]:
                        self.max_counts_percentiles[index][v][i] = self.max_counts_percentiles[index+1][v][i]
                        cont = cont + 1
        print("cont: "+str(cont))

        azzerati=0
        for index in range(len(self.contatori) - 1, -1, -1):
            for v in self.genes:
                remains = self.max_counts[index][v]
                for i in range(len(self.thresholds) - 1):
                    if remains > 0:
                        if remains > self.max_counts_percentiles[index][v][i]:
                            remains -= self.max_counts_percentiles[index][v][i]
                        else:
                            self.max_counts_percentiles[index][v][i] = remains
                            remains = 0
                    else:
                        self.max_counts_percentiles[index][v][i] = 0
                        azzerati += 1
        print("azzerati: "+str(azzerati))

        import pickle
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k) + "_max_counts"
        f = open(filename, "wb")
        pickle.dump(self.max_counts, f)
        pickle.dump(self.max_counts_percentiles, f)
        f.close()
    else:
        #Devo aggiornare ordered_percentiles
        cont=0
        for v in self.genes:
            for index in range(len(self.contatori) - 2, -1, -1):
                for i in range(len(self.thresholds) - 1):
                    a = self.ordered_percentiles[index][v][i]
                    my_values = self.ordered_percentiles[index+1][v][i]
                    if len(my_values)>0:#continua, altrimenti non c'è niente da aggiungere
                        idx = np.searchsorted(a, my_values)
                        num = np.count_nonzero(idx >= len(a))
                        if num<len(my_values):#allora almeno un valore va inserito
                            cont=cont+1
                            b = np.insert(a, idx, my_values)
                            b = b[0:self.max_counts_percentiles[index][v][i]]
                            #if (len(b) > 0):#dato che num<len(my_values) quest'altra condizione non servirebbe ma giusto per essere sicuri
                            self.ordered_percentiles[index][v][i] = b
        print("cont: " + str(cont))

        #Dai ordered_percentiles creo best_vectors
        for index in range(len(self.contatori)):
            for v in self.genes:
                b_v = np.asarray([])
                for i in range(len(self.thresholds) - 1):
                    # if len(b_v)==self.max_counts[index][v]:
                    # break
                    if self.max_counts_percentiles[index][v][i] > 0:
                        #[:self.max_counts_percentiles[index][v][i] è ripetitivo ma giusto per sicurezza (questa parte non è expensive comunque)
                        b_v = np.concatenate((b_v, self.ordered_percentiles[index][v][i][
                                                   0:self.max_counts_percentiles[index][v][i]]))
                self.best_vectors[index][v][self.k-1] = b_v
        del self.ordered_percentiles


        import pickle
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k)
        f=open(filename,"wb")
        pickle.dump(self.best_vectors,f)
        f.close()

def update_creaBestVectorsDistanza_iterations_percentiles(self):
    index=self.index
    for i in range(len(self.contatori)):
        if self.cont>=self.contatori[i]:
            index=i
    if index!=self.index:
        self.index = index
        print(self.index)


#####################################################################
###########BEST VECTORS DISTANZA 1 ITERATIONS########################
#####################################################################

def pre_creaBestVectorsDistanza_iterations(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)
    self.best_vectors = [ [ [[] for n in range(self.k)] for j in range(self.max_node + 1)]  for i in range(len(self.contatori))]
    creaBestVectors1_iterations(self)

    import pickle
    filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza1"
    f=open(filename,"wb")
    pickle.dump(self.best_vectors,f)
    f.close()
    raise ValueError

def creaBestVectors1_iterations(self):
    M=self.G.copy()
    for v in self.G.nodes:
        self.best_vectors[0][v][0]=updateBestVector_iterations(self,self.G,v)

    for index in range(1,len(self.contatori)):
        neighbors=set()
        for v in self.sorted_vertices[self.contatori[index-1]:self.contatori[index]] :
            neighbors.update(self.G.neighbors(v))
        self.G.remove_nodes_from(self.sorted_vertices[self.contatori[index-1]:self.contatori[index]])
        for v in self.G:
            if v in neighbors:
                self.best_vectors[index][v][0]=updateBestVector_iterations(self,self.G, v)
            else:
                self.best_vectors[index][v][0]=self.best_vectors[index-1][v][0]
    self.G=M

def updateBestVector_iterations(self,G,v):
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

######################################################################
###########BEST VECTORS DISTANZA 2 ITERATIONS#########################
######################################################################
def pre_creaBestVectorsDistanza_iterations(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    import pickle
    filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k - 1)
    f=open(filename,"rb")
    self.best_vectors = pickle.load(f)
    f.close()

    if self.onlyCount:
        self.max_counts =[ [0 for i in range(self.max_node+1)]  for j in range(len(self.contatori))]
    else:
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k) + "_max_counts"
        f=open(filename,"rb")
        self.max_counts=pickle.load(f)
        f.close()

    self.index=0


def crea_creaBestVectorsDistanza_iterations(self,C, vec):
    if self.onlyCount==False:
        diff = which_diff(vec)
        neighbors = set()
        for c in C:
            neighbors.update(self.G.neighbors(c))
        for u in neighbors:
            b = self.best_vectors[self.index][u][self.k-1]
            ord = np.sort(diff)
            ord = ord[0:self.max_counts[self.index][u]]
            idx = np.searchsorted(b, ord)
            b = np.insert(b, idx, ord)[0:self.max_counts[self.index][u]]
            self.best_vectors[self.index][u][self.k-1] = b
    else:
        diff = which_diff(vec)
        max = len(diff)
        neighbors = set()
        for c in C:
            neighbors.update(self.G.neighbors(c))
        for u in neighbors:
            if self.max_counts[self.index][u] < max:
                self.max_counts[self.index][u] = max


def save_creaBestVectorsDistanza_iterations(self):
    if self.onlyCount:
        cont = 0
        for index in range(len(self.contatori) - 2, -1, -1):
            for v in self.genes:
                if self.max_counts[index][v] < self.max_counts[index + 1][v]:
                    self.max_counts[index][v] = self.max_counts[index + 1][v]
                    cont=cont+1
        print("cont: "+str(cont))

        import pickle
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k) + "_max_counts"
        f = open(filename, "wb")
        pickle.dump(self.max_counts, f)
        f.close()
    else:
        #cont=0
        for v in self.genes:
            for index in range(len(self.contatori) - 2, -1, -1):
                a=self.best_vectors[index][v][self.k-1]
                my_values=self.best_vectors[index+1][v][self.k-1]
                idx=np.searchsorted(a, my_values)
                num = np.count_nonzero( idx >= len(a))
                #if num >0:
                #    cont=cont+1
                b=np.insert(a, idx, my_values)
                b=b[0:self.max_counts[index][v]]
                if(len(b)>0):
                    self.best_vectors[index][v][self.k-1]=b
        #print("cont: " + str(cont))

        import pickle
        filename = self.parameters["bestVectors"] + "/" + "BestVectorsDistanza" + str(self.k)
        f=open(filename,"wb")
        pickle.dump(self.best_vectors,f)
        f.close()

def update_creaBestVectorsDistanza_iterations(self):
    for i in range(len(self.contatori)):
        if self.cont>=self.contatori[i]:
            self.index=i

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

def pre_bound_order_improved_iterations_percentiles(self):
    self.matrix = toMatrix(self, self.G.nodes)
    ordinamentoVertici_bound_order_improved(self)

    import pickle
    filename=self.parameters["bestVectors"]
    f=open(filename ,"rb")
    self.best_vectors = pickle.load(f)
    f.close()

    """
    self.methods=["normal", "iterations", "iterations_singolo","iterations_superiore","iterations_superiore"]
    diz={"normal":0, "iterations":1, "iterations_singolo":2, "iterations_superiore":3}
    self.max_dist=len(methods)
    for i in range(self.max_dist):
        self.methods[i]=diz[self.methods[i]]
    """

    """
    if dist==3:
        best_vector = self.best_vectors[0][self.cont][dist - 1]
        #best_vector=self.best_vectors[self.index][0][dist-1]
    else:
        best_vector=self.best_vectors[self.index][v][dist-1]#già ordinato in ordine crescente
    """

def bound_order_improved_iterations_percentiles(self,C,vecC):
    dist=self.k-len(C)
    lista=which_diff(vecC)
    dec=np.asarray(lista)
    bests=[]
    if(dist<=7):
        v=C[len(C)-1]
        if dist<=2:
            best_vector = self.best_vectors[self.index][v][dist - 1]#iterations e iterations_percentiles
        else:
            if dist<=3:
                best_vector = self.best_vectors[self.index][0][dist - 1]#iterations_percentiles_singolo
            else:
                best_vector = self.best_vectors[0][self.cont][dist - 1]#Distanze superiori

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


def update_bound_order_improved_iterations_percentiles(self):
    index=self.index
    for i in range(len(self.contatori)):
        if self.cont>=self.contatori[i]:
            index=i
    if index!=self.index:
        self.index = index
        print(self.index)
    print(self.cont)
    print(self.levels)

