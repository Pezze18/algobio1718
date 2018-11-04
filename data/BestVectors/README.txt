Base: Contiene best_vectors di base(senza iterazioni) e divisi per nodi per le distanze 1 e 2 i.e. best_vectors[v][dist]
Iterations: best_vectors creato con le iterazioni(senza percentili) per distanze 1 e 2  best_vectors[index][v][dist]
IterationsPercentilesDistanza7:  best_vectors creato con le iterazioni con percentili per distanze 1 e 2, per iterazioni percentili singolo per distanza 3 e per distanza>=4 con distanze superiori
Se distanza 1 e 2 best_vectors[index][v][dist], per distanza 3 best_vectors[index][0][dist], per distanze >=4 best_vectors[0][cont][dist]
IterationsPercentilesSingoloDistanza4:  best_vectors creato con le iterazioni con percentili per distanze 1 e 2, per iterazioni percentili singolo per distanza 3 e 4 e per distanza>4 con distanze superiori

parametri per eseguire bound_order 
crea e onlyCount vanno messi false mentre il parametro si imposta il file BestVectors che si vuole eseguire mentre il metodo è bound_order_improved_iterations_percentiles
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
    'bound': True,          #creaBestVectorsDistanza1_iterations_percentiles
    'method': "bound_order_improved_iterations_percentiles",
    'bestVectors':"../data/BestVectors/BestVectorsDistanza7",
    "crea":False,
    "onlyCount":False
}

parametri per creare vettori
Durante la creazione si imposta il parametro bestVectors in cui si imposta il path di dove salvare il BestVector, per esempio "../data/BestVectors/"
Per k=1: si invoca con k=1 creaBestVectorsDistanza1_iterations_percentiles
Per k=2: si imposta k=2 e si invoca creaBestVectorsDistanza_iterations_percentiles e si imposta crea=True e onlyCount=True e crea il BestVectorsK_max_counts
          dopodichè si imposta onlyCount=False e si riesegue creando il BestVectorsK
Per k=3 e k=4 si invoca creaBestVectorsDistanza_iterations_percentiles_singolo
Per k>4 si invoca creaBestVectorsDistanzeSuperiori
  
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
    'bound': True,          #creaBestVectorsDistanza1_iterations_percentiles
    'method': "bound_order_improved_iterations_percentiles",
    'bestVectors':"../data/BestVectors/BestVectorsDistanza7",
    "crea":False,
    "onlyCount":False
}