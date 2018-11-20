import pandas as pd
import os
import warnings



#########WARNING##################################################
# Tutte le esecuzioni devono essere riuscite(no failure)
# possibilmente lo stesso valore k per un metodo ci sia una volta sola
######################################################################

strategy = "enumerate"
transform = True
ordina = False
if(strategy=="enumerate"):
    names_order = [
        "oldVersion",
        "noboundMatrix",
        "noboundLevelsVec",
        "nobound_migliorato",
        "bound_min",
        "bound_min_migliorato",
        "bound_min_migliorato_optimized",
        "bound_min_migliorato_iterations",
        "bound_min_migliorato_iterations_optimized",
        "bound_order",
        "bound_order_optimized",
        "bound_order_improved_iterations_percentiles",
        "bound_means",
        "bound_means_optimized",
        "bound_means_iterations",
        "bound_means_iterations_optimized"
    ]
else:
    names_order=[
                 "oldVersion",
                 "onlyBFS",
                 "onlyBFSAndProbCover",
                 #"onlyBFSAndComplement",
                 "BFSAndProbCoverAndComplement",
                 "BFSAndLevelsVec"]


def prepareNumber(value):
    # print("value: "+str(value))
    final = ""
    length = len(value)
    number = int(length / 3.0)
    part = length - 3 * number
    # print(part)

    if (part > 0):
        final = final + value[0:part]
        if (number >= 1):
            final = final + "." + value[part:part + 3]
    else:
        if (
                number >= 1):  # if superfluo, se part==0 allora deve esserci almeno una terna, altrimenti numero sarebbe vuoto assurdo !
            final = final + value[0:3]
    for i in range(1, number):
        j = i * 3 + part
        final = final + "." + value[j] + value[j + 1] + value[j + 2]
    # print("value_final: "+str(final))
    return final


def prepareNumbers(lista):
    # print("Before: "+str(lista))
    lista = lista[1:len(lista) - 1]
    lista = lista.replace(" ", "")
    lista = lista.split(",")
    for i in range(len(lista)):
        # print(str(i)+":"+lista[i])
        lista[i] = prepareNumber(lista[i])
    # print(lista)
    # print("Before2: "+str(lista))
    final = "["
    for i in range(len(lista)):
        final += lista[i] + ", "
    final = final[0:len(final) - 2]
    final += "]"
    # print("After: "+str(final))
    return final


def convertTime(time):
    (time, seconds) = time.split("m")
    time = int(time)
    seconds = (float(seconds.split("s")[0]))

    ms = seconds - int(seconds)
    seconds = int(seconds)

    h = int(time / 60)
    m = time - 60 * h
    if h == 0:
        h = "00"
    if m == 0:
        m = "00"
    if seconds == 0:
        seconds = "00"

    h = str(h)
    m = str(m)
    s = str(seconds)
    ms = str(ms * 1000)[0:3]

    if len(h) == 1:
        h = "0" + h
    if len(m) == 1:
        m = "0" + m
    if len(s) == 1:
        s = "0" + s

    # return str(h)+"h:"+str(m)+"m:"+str(seconds)+"s:"+str(ms)+"ms"
    return str(h) + ":" + str(m) + ":" + str(seconds)


def sortSolution(sol):
    sol = sol[1:len(sol) - 1]
    sol = sol.replace("'", "").replace("[", "").replace("]", "")
    sol = sol.split(",")
    #print("Before:")
    #print(sol)
    sol.sort()
    #print("After:")
    #print(sol)
    #print("________________")
    return str(sol).replace("'", "")






columns = ["method"] + ["k=" + str(i) for i in range(1, 11)]
dfTimes = pd.DataFrame([], columns)
rows_list_times = []
rows_list_levels = []
rows_list_scores = []
rows_list_solutions = []

for method in os.listdir("out/" + strategy):
    print(method)
    rigaTimes = ["" for k in range(0, 11)]
    rigaTimes[0] = method
    rigaLevels = ["" for k in range(0, 11)]
    rigaLevels[0] = method
    rigaScores = ["" for k in range(0, 11)]
    rigaScores[0] = method
    rigaSolutions = ["" for k in range(0, 11)]
    rigaSolutions[0] = method

    for directory in os.listdir("out/" + strategy + "/" + method):
        #print(directory)
        if directory.startswith("k="):
            end = directory.find("_")
            if (end == -1):
                end = len(directory)
            k = int(directory[2:end])
            # print("k: "+str(k))
            for file in os.listdir("out/" + strategy + "/" + method + "/" + directory):
                path = "out/" + strategy + "/" + method + "/" + directory + "/" + file
                # Manage row
                if file.startswith("k="):
                    if os.path.getsize(path) == 0:
                        warnings.warn(path + " is empty")
                    else:
                        f = open(path)
                        lines = []
                        ultimo_index = -1
                        index = 0
                        for line in f:
                            lines.append(line)
                            if line.find("_________________") != -1:
                                ultimo_index = index
                            index += 1
                        # print(ultimo_index)
                        #print(path)
                        #print(lines)

                        FinalsolutionID = lines[ultimo_index - 5]
                        FinalsolutionNames = lines[ultimo_index - 4]
                        FinalsolutionScore = lines[ultimo_index - 2]
                        ElapsedTime = lines[ultimo_index - 1]
                        Levels = lines[ultimo_index - 7]
                        Levels = "[" + Levels[4:len(Levels)].replace("\n", "")

                        # print(FinalsolutionID)
                        # print(FinalsolutionNames)
                        # print(FinalsolutionScore)
                        # print(ElapsedTime)
                        # print(Levels)

                        if (strategy == "enumerate"):
                            if (transform):
                                Levels = prepareNumbers(Levels)
                            rigaLevels[k] = Levels
                        # rigaGenes[k]=FinalsolutionNames.split(":")[1].replace("\n","")
                        # rigaGenes[k]=rigaGenes[k][2:len(rigaGenes[k])]
                        # print(rigaGenes)

                        rigaScores[k] = FinalsolutionScore.split(":")[1].replace(" ", "")

                        FinalsolutionNames = FinalsolutionNames.split(":")[1].replace(" ", "")
                        if (transform):
                            FinalsolutionNames = sortSolution(FinalsolutionNames)
                        rigaSolutions[k] = FinalsolutionNames

                        f.close()

                if file.startswith("commands.job.e"):
                    f = open(path)
                    for line in f.readlines():
                        if line.startswith("real"):
                            time = line[4:len(line)]
                            time = time.strip()
                            # print("time:"+str(time))
                            if (transform):
                                time = convertTime(time)
                            rigaTimes[k] = time
                    f.close()

    rows_list_times.append(rigaTimes)
    rows_list_levels.append(rigaLevels)
    rows_list_scores.append(rigaScores)
    rows_list_solutions.append(rigaSolutions)

# ordina
if (ordina):
    lista = [rows_list_times, rows_list_levels, rows_list_scores, rows_list_solutions]
    for j in range(len(lista)):
        rows = lista[j]
        names = []
        indexes_order = []
        # print(rows)
        for row in rows:
            # print(row)
            names.append(row[0])
        # print(names)
        for name in names_order:
            #print(name)
            indexes_order.append(names.index(name))

        lista_rows = []
        for i in indexes_order:
            lista_rows.append(lista[j][i])
        lista[j] = lista_rows

    rows_list_times = lista[0]
    rows_list_levels = lista[1]
    rows_list_scores = lista[2]
    rows_list_solutions = lista[3]

dst = "data/summary_tables/"
dfTimes = pd.DataFrame(rows_list_times)
dfTimes.columns = columns
dfTimes.to_csv(dst + "times_" + strategy + ".csv")

if (strategy == "enumerate"):
    dfLevels = pd.DataFrame(rows_list_levels)  # ,columns
    dfLevels.columns = columns
    dfLevels.to_csv(dst + "levels_" + strategy + ".csv")

dfScores = pd.DataFrame(rows_list_scores)
dfScores.columns = columns
dfScores.to_csv(dst + "scores_" + strategy + ".csv")

dfSolutions = pd.DataFrame(rows_list_solutions)
dfSolutions.columns = columns
dfSolutions.to_csv(dst + "solutions_" + strategy + ".csv")

