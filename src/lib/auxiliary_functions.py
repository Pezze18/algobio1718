import numpy as np
import math

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

def prob_cover_max_version(self,C):#numpy #min version
    v=vectorization_solution(self,C)
    som = np.sum(v)
    return len(self.samples)-som

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
def BFS_complete_node(self,radice,creaLevels=True,creaPredecessori=True,creaVec=True,creaEtichetta=True, creaDepths=True):
    visit = [False for i in range(10000)]
    visit[radice]=True

    #if creaLevels:
    L=[[] for i in range(self.k+1)]#0,1,...k
    L[1]=[radice]#al livello 1 c'è la radice e k è l'ultimo livello
    self.L=L

    if creaPredecessori:
        pred = [0 for i in range(10000)]# se alla fine pred[g] è ancora a 0, allora g non è raggiungibile da v
        pred[radice]=radice
        self.pred= pred
    if creaVec:
        shortestVec = [0 for i in range(10000)]
        shortestVec[radice] = np.ones(len(self.samples))
        self.shortestVec=shortestVec
    if creaEtichetta:
        labels=[False for i in range(10000)]
        labels[radice]=True# v viene sicuramente scelta in C_v
        self.labels=labels
    if creaDepths:
        depths=[-1 for i in range(10000)]
        depths[radice]=0
        self.depths=depths

    for k in range(2,self.k+1):
        #print("Livello: "+str(k))
        for g in L[k-1]:
            neighbors=self.G.neighbors(g)
            for u in neighbors:
                if visit[u] == False:
                    visit[u] = True
                    if creaDepths:
                        depths[u] = k
                    if creaLevels:
                        L[k].append(u)
                    if creaPredecessori:
                        self.pred[u] = g
                    if creaVec:
                        shortestVec[u] = np.multiply(shortestVec[g], self.matrix[u])
        #if creaLevels==False:
        #    L[k-1]=0
    #if creaLevels == False:
    #    del self.L

def findAncestor(self,v,s):#newC è l'insieme complemento e father è l'ancestor comune(nel caso peggiore == radice)
    newC=[s]
    #print("v: "+str(v))
    #print("s: " + str(s))
    #print("father: "+str(self.pred[s]))
    while True:
        father=self.pred[s]
        #print("father: " + str(father))
        #if father==v:
        #    return newC,v
        if self.labels[father]:#allora è già stato selezionato
            #print("Trovato ancestor")
            return newC,father#Nel caso peggiore ritorna la radice(primo nodo selezionato)
        s=father
        newC.append(s)


def BFS_complete(self,creaLevels=True,creaPredecessori=True,creaVec=True,creaEtichette=True,creaDepths=True):
    if creaLevels:
        L=[ [[] for v in range(self.k)] for i in range(10000)]
    visit = [[False for v in range(10000)] for i in range(10000)]
    if creaPredecessori:
        pred = [[None for v in range(10000)] for i in range(10000)]
    if creaVec:
        shortestVec = [[None for v in range(10000)] for i in range(10000)]
    if creaEtichette:
        labels=[[False for v in range(10000)] for i in range(10000)]#100000*100000= 10 000 000 000

    for v in self.G:
        L[v][0]=v
        pred[v][v]=v
        shortestVec[v][v]=self.matrix[v]#per il percorso v, il vettore per v(radice) è self.matrix[v]
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

