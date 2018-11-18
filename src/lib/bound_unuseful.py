import networkx as nx
from lib.core import *
import numpy as np
import math
#import tensorflow as tf





###################################
####### BOUND_PROD ########
###################################
def pre_bound_prod(self):
    self.matrix = toMatrix(self.samples, self.G.nodes)

    #cerco min
    min=10
    for g in self.matrix:
        for n in self.matrix[g]:
            if(n<min and n!=0):
                min=n
    print("min: "+str(min))

    #aggiorni gli 0 con min
    for g in self.matrix:
        for i in range(0,len( self.matrix[g])):
            if(self.matrix[g][i]==0):
                self.matrix[g][i]=min


def bound_prod(self,C):#DA COMPLETARE
    return False

def ordinamentoVertici_bound_prod(self):
    sorted_vertices = [t[0] for t in sorted(self.G.degree, key=lambda t: t[1])]
    return sorted_vertices

def update_bound_prod(self):
    return False


###################################
####### BOUND_ALTERNATIVE ########
###################################
""""
debugging in inspectNode function quando len(C)==k
 score=1/score
 print("ActualScore: "+str(score))
 print("BoundScore: " + str(prob_cover_alternative(self,C)))
 print("rapporto: "+str())
 print()
 """
def pre_bound_alternative(self):
    self.matrix = toMatrix(self.samples, self.G.nodes)

    self.alternative_matrix={}
    for g in self.matrix:
        self.alternative_matrix[g]=np.exp(-self.matrix[g])

def bound_alternative(self,C):
    dist=self.k-len(C)
    return False

def ordinamentoVertici_bound_alternative(self):
    sorted_vertices = [t[0] for t in sorted(self.G.degree, key=lambda t: t[1])]
    return sorted_vertices

def update_bound_alternative(self):
    return False

###################################
####### BOUND_LOG ########
###################################
def pre_bound_log(self):
    self.matrix = toMatrix(self.samples, self.G.nodes)

    #cerco min
    min=10
    for g in self.matrix:
        for n in self.matrix[g]:
            if(n<min and n!=0):
                min=n
    print("min: "+str(min))

    #aggiorni gli 0 con min
    for g in self.matrix:
        for i in range(0,len( self.matrix[g])):
            if(self.matrix[g][i]==0):
                self.matrix[g][i]=min

    #Scalatura

    t=1#1/min
    for g in self.matrix:
        for i in range(0,len( self.matrix[g])):
            self.matrix[g][i]=t*self.matrix[g][i]

    self.sum={}
    for g in self.matrix:
        #print(np.log(self.matrix[g]))
        self.sum[g]=np.max(np.log(self.matrix[g]))
    print(self.sum)


def bound_log(self,C):
    return False
    dist=self.k-len(C)
    if(dist==1):
        for u in self.G.neighbors(C[0]):
            bestS=self.sum[C[0]]+self.sum[u]
            e=2.71828
            print("BestS: "+str(bestS))
            bestS=math.exp(bestS)
            print("BestS: " + str(bestS))
            print()
            if(bestS>self.best_score):
                print("Hello")
                return True
    return False

def ordinamentoVertici_bound_log(self):
    sorted_vertices = [t[0] for t in sorted(self.G.degree, key=lambda t: t[1])]
    return sorted_vertices

def update_bound_log(self):
    return False


###################################
####### BOUND_KANTOROVICH ########
###################################
"""
debuggin per kantorovich nella funzione inspectNode con len(C)==k
x=vectorization_solution( self, C[0:len(C)-1])
y=vectorization_solution( self, [C[len(C)-1]])
rec=np.reciprocal(y)
rap=np.multiply(x,rec)
M=np.max(rap)
m=np.min(rap)
A=(m+M)/2
G=math.sqrt(m*M)
print("m: "+str(m))
print("M: " + str(M))
optimalScore=np.linalg.norm(x)*np.linalg.norm(y)*(G/A)
score = self.scoring_function(self,C)
print("OptimalScore: "+str(optimalScore))
print("ActualScore: "+str(score))
print()
"""

def pre_bound_kantorovich(self):
    print("Inizio Ordinamento")
    self.matrix = toMatrix(self.samples, self.G.nodes)

    self.best_subsets=[ ["PTEN"], ["PTEN", "TP53"],["PTEN", "TP53", "RB1CC1"], ["PTEN", "TP53", "PTK2","EGFR"] ]
    self.best_subsets=[[self.str_to_id[x] for x in sol ] for sol in self.best_subsets]
    print(self.best_subsets)
    self.best_scores=[prob_cover(self,sol) for sol in self.best_subsets]
    print(self.best_scores)
    self.best_scores=[[-1]]+self.best_scores
    self.best_subsets=[[]]+self.best_subsets

    rapporto={}
    for g in self.matrix:
        rapporto[g]=np.min(self.matrix[g])/np.max(self.matrix[g])
        #print(str(rapporto[g]) +" "+str(np.min(self.matrix[g])))

    self.best_rapporti={}

    self.visit=[False for i in range(0,len(self.G.nodes)+1)]

    for g in self.G:
        L=BFS(self, g)
        #print(L)
        self.best_rapporti[g]={}
        for k in range(1,self.k):
            lista=[]
            for j in range(1,k):
                lista=lista+L[j]
            lista=[rapporto[x] for x in lista  ]
            #print(lista)
            lista.sort()
            lista=lista[0:k]

            self.best_rapporti[g][k]=np.prod(lista)
    """"
    for g in self.G:
        for k in range(1,self.k):
            self.best_rapporti[g][k]*= ( math.sqrt(self.best_scores[k])/math.sqrt(math.sqrt(len(self.samples)))  )
    del self.visit
    """
    print("Fine Ordinamento")

def bound_kantorovich(self,C):
    vec=vectorization_solution(self,C)
    norma=np.linalg.norm(vec,None)
    dist=self.k-len(C)
    min=np.min(vec)
    max=np.max(vec)
    for u in C:
        #print("norma di x: "+str(norma))
        rapporti=math.sqrt(self.best_scores[dist])/math.sqrt(math.sqrt(len(self.samples)))
        #print("bound norma y: "+str(rapporti))
        #print("best rapporto y: "+str(self.best_rapporti[u][dist]))
        #print("rapporto di x: "+str(min/max))
        bestS=self.best_rapporti[u][dist]*norma*(min/max)*rapporti
        #print("bestS: "+str(bestS))
        if(bestS>self.best_score):
            return True

    return False


def ordinamentoVertici_bound_kantorovich(self):
    sorted_vertices = [t[0] for t in sorted(self.G.degree, key=lambda t: t[1])]
    return sorted_vertices

def update_bound_kantorovich(self):
    return True