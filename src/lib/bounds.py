import numpy as np


###################################
########### NO_BOUND #############
###################################
def pre_nobound(self):
    self.matrix = toMatrix(self, self.G.nodes)
    #Ordinamento Vertici
    self.sorted_vertices= [t[0] for t in sorted(self.G.degree, key=lambda t: t[1])]

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
        if len(samples[p].intersect(C))>0
    ])

def set_cover(self,C):
    samples = self.samples
    return set([
        p for p in samples
        if len(samples[p].intersect(C)) > 0
    ])