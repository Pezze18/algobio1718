import networkx as nx
import lib.bounds as bounds
from lib.bounds import *
import numpy as np
import math
import bottleneck as bottle
import lib.combinatorial as combinatorial
from lib.combinatorial import *

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
        self.levelsVecUse=True#parameters["levelsVec"]
        self.max_node=9859
        self.crea=parameters["crea"]
        self.onlyCount=parameters["onlyCount"]

        self.pre=getattr(bounds,"pre_"+parameters["method"])
        self.pre(self)

        if self.prob and self.bound:
            self.bound_function = getattr(bounds, parameters["method"])
            self.update_function = getattr(bounds, "update_" + parameters["method"])
        if self.crea:
            self.crea_function = getattr(bounds, "crea_" + parameters["method"])
            self.save_function = getattr(bounds, "save_" + parameters["method"])

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

            if self.crea:
                self.crea_function(self, C, self.levelsVec[size])

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

        if self.crea:
            self.save_function(self)

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
        self.genes=list(G.nodes)
        self.parameters=parameters

        if(self.prob):
            self.matrix = toMatrix(self, self.G.nodes)
        if (self.prob and self.bound == False):
            self.comb_alg=getattr(combinatorial, "combinatorial_algorithm_prob_"+parameters["method"])
        if (self.prob and self.bound):
            self.bound_function=getattr(bounds,parameters["method"])
            self.pre_function=getattr(bounds,"pre_"+parameters["method"])

    def combinatorial_algorithm(self):
        if(self.prob and self.bound==False):
            return self.comb_alg(self)
        if (self.prob and self.bound):
            self.pre_function(self)
            return combinatorial_algorithm_prob_BFSAndLevelsVecAndPruning(self)
        return combinatorial_algorithm_det(self)


class FloydWarshall:
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
    def execute(self):
        d=[[10000 for j in range(10000)] for i in range(10000)]#costo infinito per tutti inizialmente
        for i in self.G.nodes:
            for j in self.G.nodes:
                if i!=j:
                    d[i][j]=np.multiply(self.matrix[i], self.matrix[j])

        nodes=set(self.G.nodes)
        for h in range(10000):
            if h not in nodes:
                continue
            for i in nodes:
                for j in nodes:
                    d_ihhj=np.multiply(d[i][h],d[h][j])
                    JoinScore=np.sum(d_ihhj)
                    if JoinScore < d[i][j]:
                        d[i][j]=d_ihhj



class UpperBoundEstimate:
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
    def execute(self):
        percentiles=[i*10 for i in range(0,10)]
        totale=[]
        for g in self.G.nodes:
            totale+=list(self.matrix[g][self.matrix[g]<1])
        thresholds=np.percentile(totale,percentiles)
        thresholds=list(thresholds)+[1]
        print(thresholds)

        orderedMatrix = [[] for i in range(10000)]
        max_counts = [0 for i in range(10000)]  # mi pare sia corretto
        for g in self.G.nodes:
            orderedMatrix[g] = list(which_diff(self.matrix[g]))
            max_counts[g]= len(orderedMatrix[g]) # ad ogni nodo è associato il max_count




        max_counts_percentiles = [0 for i in range(10)]
        for g in self.G.nodes:
            counts = list(np.histogram(orderedMatrix[g], thresholds)[0])
            for i in range(len(thresholds)-1):
                if max_counts_percentiles[i]<counts[i]:
                    max_counts_percentiles[i]=counts[i]
        print(max_counts_percentiles)

        ordered_percentiles=[[] for i in range(len(thresholds)-1)]
        cont=0
        for g in self.G.nodes:
            indices=list(np.digitize(orderedMatrix[g],thresholds)-1)
            cont+=1
            for i in range(len(indices)):
                ordered_percentiles[indices[i]].append(orderedMatrix[g][i])
        #print("cont: "+str(cont))
        for i in range(len(ordered_percentiles)):
            ordered_percentiles[i] = bottle.partition(ordered_percentiles[i], max_counts_percentiles[i])[:max_counts_percentiles[i]]

        #print(ordered_percentiles)

        counts_k = sorted(max_counts,reverse=True)[0:self.k]
        best_vectors=[[] for i in range(self.k)]
        for i in range(len(counts_k)):
            index_perc = len(thresholds) - 2
            best_vectors[i]=np.ones(len(self.samples)-counts_k[i])
            remains=counts_k[i]
            while(remains>0):
                if remains>len(ordered_percentiles[index_perc]):

                    best_vectors[i]=np.concatenate((best_vectors[i],ordered_percentiles[index_perc]))
                    remains-=len(ordered_percentiles[index_perc])
                else:
                    best_vectors[i]=np.concatenate((best_vectors[i],
                                                   ordered_percentiles[index_perc][0:remains]))
                    remains=0
        result=best_vectors[0]
        for i in range(1,len(best_vectors)):
            result=np.multiply(result,best_vectors[i])

        score=np.sum(result)
        score_max=len(self.samples)-score
        print("MinVersion: "+str(score))
        print("MaxVersion: "+str(score_max))



