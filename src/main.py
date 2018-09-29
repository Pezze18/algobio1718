from lib.inout import *
from lib.core import *
from lib.bounds import *
import time


def main():
    parameters = parse_command_line()
    k = parameters['k']
    delta = parameters['delta']
    strategy = parameters['strategy']
    G = load_network(parameters['genes_input'])
    id_to_str, str_to_id = make_maps(parameters['proteins_input'])
    parameters["str_to_id"]=str_to_id
    parameters["id_to_str"]=id_to_str
    filter = read_filter(parameters['filter_input'])

    t_start = time.time()
    if parameters['prob']:
        patients = read_patients_prob(parameters['samples_input'], str_to_id, filter=filter)
    else:
        patients = read_patients(parameters['samples_input'], str_to_id, filter=filter)

    if strategy == 'combinatorial':
        if parameters['prob']:
            C, score_C = prob_combinatorial_algorithm(G,k,patients)
        else:
            C, P_C = combinatorial_algorithm(G,k,patients)
            score_C = len(P_C)
    elif strategy == 'enumerate':
        BDDE_obj = BDDE(G, patients, k, parameters)
        BDDE_obj.enumeration_algorithm()

        C = BDDE_obj.best_subgraph
        score_C = BDDE_obj.best_score
    else:
        raise ValueError("Unkown strategy given. Input: " + str(strategy))

    t_end = time.time()
    print("_________________")
    print("Final solution (ids): ", C)
    print("Final solution (genes' names): ", [id_to_str[el] for el in C])
    if parameters['prob']:
        print("Final solution score: ", score_C)
        print("Final solution score(max version): ", len(patients)-score_C)
    else:
        print("Final solution cardinality: ", score_C)
    print("Elapsed time: ", time.strftime("%H:%M:%S", time.gmtime(t_end - t_start)))
    print("_________________")

if __name__ == "__main__":
    main()