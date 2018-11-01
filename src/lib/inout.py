import pandas as pd
import argparse
import os
import time
import datetime
import networkx as nx

def check_parameters(parameters):
    #Important check
    if (parameters['prob'] == False and parameters['bound']):
        raise ValueError("you cannot execute bound option without abilitate prob option")
    if (parameters['prob'] and parameters['method']=="det"):
        raise ValueError("you cannot execute prob mode using method det")
    if (parameters["bound"] and parameters['method'] == "nobound"):
        raise Exception("You cant use bound option when you execute nobound")

    ##Check files path##
    if not os.path.isfile(parameters['filter_input']):
        #for file in os.listdir("../../../../"):
        #    print(file)
        raise FileNotFoundError("Can't find the 'proteins' file; filename given: " + parameters['filter_input'])

    if not os.path.isfile(parameters['proteins_input']):
        raise FileNotFoundError("Can't find the 'proteins' file; filename given: " + parameters['proteins_input'])

    if not os.path.isfile(parameters['samples_input']):
        raise FileNotFoundError("Can't find the 'samples' file; filename given: " + parameters['samples_input'])

    if not os.path.isfile(parameters['genes_input']):
        raise FileNotFoundError("Can't find the 'genes' file; filename given: " + parameters['genes_input'])

    #Other checks
    if not (0 < parameters['delta'] <= 1):
        raise ValueError("The given 'delta' is not a valid  value (i.e. in ]0,1]). Given: " + str(parameters['delta']))

    k=parameters['k']
    if not (k >= 1) or type(k) != int:
        raise ValueError("The given 'k' is not a valid value (i.e. integer greater than 1). Given: " + str(k))

# Parameters of default
matriceProb="matriceProb.csv"
matriceBinaria="matriceBinaria.csv"
parameters = {
    'k':5,
    'proteins_input': "../data/hint+hi2012_index_file.txt",
    'samples_input': "../data/",
    'genes_input': "../data/hint+hi2012_edge_file.txt",
    'filter_input': '../data/mutated_expressed_genes.txt',
    'delta': 0.8,
    'prob': True,
    'strategy': 'enumerate',  # options: enumerate,combinatorial
    'best_score': 10000000,  # maximum for a single gene
    'bound': True,
    'method': "bound_order_improved_iterations",#bound_order_improved#creaBestVectorsDistanza2
    'bestVectors':"../data/bestVectors2",
    "crea":False,
    "onlyCount":False
}
# options method enumerate: det, nobound, bound_min, bound_min_migliorato, bound_min_migliorato_iterations
# bound_fast, bound_kantorovich, bound_means, bound_means_iterations
#options method combinatorial: onlyBFS BFSAndProbCover BFSAndLevelsVec

if parameters['prob']:
    parameters['samples_input']+=matriceProb
else:
    parameters['samples_input']+=matriceBinaria
#check_parameters(parameters)#controllo che i dati siano idonei/validi <-- spostare in main e aggiungere il parametro per local o cluster_script

def parse_command_line():
    """
    py:function:: parse_command_line()
    Parses the command line.
    :return: parameters, which is a dictionary. i.e. parameters['parameter'] = value_of_such_parameter.
    """

    #Gestisce il parser, specifica tutti i possibili campi passabili tramite linea di comando in input
    args=handleParser()
    update_parameters(args, parameters) #aggiorno il dizionario di parametri di default con campi aggiornat
    check_parameters(parameters)        # Leggo i dati ricevuti, verifico che siano idonei

    #Debuggging
    #print(parameters)
    return parameters


def handleParser():
    parser = argparse.ArgumentParser(description='Parser to receive input parameters.')

    parser.add_argument('--proteins_input', type=str, help='File path to the gene-to-id data set')
    parser.add_argument('--samples_input', type=str, help='File path to the samples data set')
    parser.add_argument('--genes_input', type=str, help='File path to the genes (graph) data set')
    parser.add_argument('--filter_input', type=str, help='File path to the filter data set')
    parser.add_argument('--delta', type=float, help='Delta parameter: edge weight tolerance')
    parser.add_argument('--k', type=int, help='Cardinality of the solution to be returned')
    parser.add_argument('--strategy', choices=['combinatorial', 'enumerate'], help='Choose the strategy')
    parser.add_argument('--prob', action="store_true", help='Type --prob if you want to use the' +
                                                            'probaiblistic version of the problem')
    parser.add_argument('--outfolder', type=str, help='The name (only!) of the folder to be created inside' +
                                                      'the \'$project_path\'/out/ directory')
    parser.add_argument('--best_score', type=float, help='The previous (k-1) best solution score')
    parser.add_argument('--bound', action="store_true", help='Type --bound if you want to use the bounds')

    pre="Method to apply for pre part(init)"
    ord=pre+", ordinamentoVertici(order list of vertices)"
    bound_f=ord+", bound_function(only if bound parameter is present)"
    update=bound_f+", update(update function to use)"
    parser.add_argument('--method', type=str, help=update )

    parser.add_argument('--scoring_function', type=str, help='Scoring function to use inside BDDE ')
    parser.add_argument('--bestVectors', type=str, help='bestVectors path for bound_order ')
    parser.add_argument('--levelsVec', action="store_true", help='Type --levelsVec if you want to use levelsVec')
    return parser.parse_args()


def update_parameters(args,parameters):
    if args.prob:
        parameters['prob'] = True

    if args.bound:
        parameters['bound'] = True

    if args.best_score:
        parameters['best_score'] = float(args.best_score)

    if args.outfolder:
        parameters['outfolder'] = args.outfolde

    if args.strategy:
        parameters['strategy'] = args.strategy

    if args.k:
        parameters['k'] = args.k

    if args.delta:
        parameters['delta'] = args.delta

    if args.genes_input:
        parameters['genes_input'] = args.genes_input

    if args.samples_input:
        parameters['samples_input'] = args.samples_input

    if args.proteins_input:
        parameters['proteins_input'] = args.proteins_input

    if args.filter_input:
        parameters['filter_input'] = args.filter_input

    if args.method:
        parameters["method"]=args.method

    if args.scoring_function:
        parameters['scoring_function'] = args.scoring_function

    if args.bestVectors:
        parameters['bestVectors']=args.bestVectors
    return parameters

def read_patients_prob(filename, genes_map, filter=set()):
    frame = pd.read_csv(filename,sep="\t")
    patients = {}
    for index,row in frame.iterrows():
        patients.setdefault(row["sampleid"], {})
        if row["geneid"] in genes_map and row["geneid"] in filter:
            patients[row["sampleid"]][genes_map[row["geneid"]]]=row["prob"]

    patients = {p:patients[p] for p in patients if len(patients[p]) > 0}

    return patients


# filtro e' un set formato dai nomi dei geni come le chiavi di genes_map quindi TTN,ecc.
def read_patients(filename, genes_map, filter=set()):

    with open(filename, "r") as file:
        patients = {}
        lines = file.readlines()
        for l in lines:
            s=l.split("\t")
            # Nota: se un gene risulta mutato in un paziente ma non compare nell'elenco dei geni dato in input,
            #       lo eliminiamo nell'elenco dei geni mutati di tale paziente, considerandolo (evidentemente)
            #       non rilevante.
            patients[s[0]]=set([genes_map[gene] for gene in s[1:] if gene in genes_map and gene in filter])
    return patients


def read_filter(filename):
    filter = set()
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            filter.add(str(line))
    return filter


def make_maps(filename):

    id_to_str={}  # dato l'id del gene, d[id] = "stringa del gene"
    str_to_id={}  # data la stringa del gene, d[stringa] = "id del gene"

    with open(filename, "r") as file:
        lines=file.readlines()
        for l in lines:
            s=l.split(" ")
            gene_id = int(s[0])
            gene_str = s[1].split("\t")[0]
            id_to_str[gene_id] = gene_str
            str_to_id[gene_str] = gene_id

    return id_to_str, str_to_id

def load_network(filename):
    with open(filename, "r") as file:
        G=nx.Graph()
        lines=file.readlines()
        for l in lines:
            s=l.split(" ")
            s=[ int(ss) for ss in s]
            G.add_edge(s[0], s[1] ,weight=s[2])

    return G

"""
    OLD/TESTING FUNCTIONS
"""

"""
def g_input_read():
    # generiamo un grafo completo un po' a caso
    G = nx.complete_graph(100)

    for (u,v) in G.edges():
        G[u][v]['weight'] = np.random.random() # ]0,1[

    return G

def p_input_read():
    p = {}
    genes = [i for i in range(100)]
    for i in range(100):
        id = "Patient#"+str(i+1)
        patient_genes = np.random.choice(genes, size=int(np.random.random()*100), replace=False)
        p[id] = set(patient_genes)

    return p


def testing(filename, genes_map):
    frame = pd.read_csv(filename,sep="\t")

    patients={}
    for index,row in frame.iterrows():
        patients.setdefault(row["sampleid"], set())
        if row["geneid"] in genes_map:
            patients[row["sampleid"]].add(genes_map[row["geneid"]])

    return patients
"""

