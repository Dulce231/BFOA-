import matplotlib.pyplot as plt
from fastaReader import fastaReader
import random
import numpy as np
import copy
from evaluadorBlosum import evaluadorBlosum

# Clase bacteria
class bacteria:
    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction = 0
        self.NFE = 0

    def showGenome(self):
        for seq in self.matrix.seqs:
            print(seq)

    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = np.array(copy.deepcopy(self.matrix.seqs))
        return newBacteria

    def tumboNado(self, numGaps):
        self.cuadra()
        matrixCopy = copy.deepcopy(self.matrix.seqs).tolist()
        gapRandomNumber = random.randint(0, numGaps)
        for _ in range(gapRandomNumber):
            seqnum = random.randint(0, len(matrixCopy) - 1)
            pos = random.randint(0, len(matrixCopy[0]))
            matrixCopy[seqnum] = matrixCopy[seqnum][:pos] + "-" + matrixCopy[seqnum][pos:]
        self.matrix.seqs = np.array(matrixCopy)
        self.cuadra()
        self.limpiaColumnas()

    def cuadra(self):
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-" * (maxLen - len(seq[i]))
        self.matrix.seqs = np.array(seq)

    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True

    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1

    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos + 1:]

    def getColumn(self, col):
        column = [self.matrix.seqs[i][col] for i in range(len(self.matrix.seqs))]
        return column

    def autoEvalua(self):
        evaluador = evaluadorBlosum()
        score = 0
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            gapCount = column.count("-")
            column = [x for x in column if x != "-"]
            pares = self.obtener_pares_unicos(column)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])
            score -= gapCount * 2
        self.blosumScore = score
        self.NFE += 1

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i + 1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)

# Clase chemiotaxis
class chemiotaxis:
    def __init__(self):
        self.parcialNFE = 0

    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 0.0
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore) ** 2.0
            total += d * np.exp(w * diff)
        return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
        return attract + repel

    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction

    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        self.parcialNFE = 0
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
            self.parcialNFE += bacteria.NFE
            bacteria.NFE = 0

    def eliminarClonar(self, path, poblacion):
        poblacion.sort(key=lambda x: x.fitness)
        for i in range(int(len(poblacion) / 2)):
            del poblacion[0]
        clones = self.clonacion(path, poblacion)
        poblacion.extend(clones)

    def clonacion(self, path, poblacion):
        poblacionClones = []
        best = max(poblacion, key=lambda x: x.fitness)
        for bacteria in poblacion:
            newBacteria = bacteria.clonar(path)
            mutacion = int((best.fitness - bacteria.fitness) / 10)
            newBacteria.tumboNado(mutacion)
            newBacteria.autoEvalua()
            poblacionClones.append(newBacteria)
        return poblacionClones

    def randomBacteria(self, path):
        bact = bacteria(path)
        bact.tumboNado(random.randint(1, 10))
        return bact

    def insertRamdomBacterias(self, path, num, poblacion):
        for _ in range(num):
            poblacion.append(self.randomBacteria(path))
            poblacion.sort(key=lambda x: x.fitness)
            del poblacion[0]

# Función clonaBest
def clonaBest(veryBest, best):
    veryBest.matrix.seqs = np.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

# Variables y parámetros
poblacion = []
path = "multiFasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)
tempBacteria = bacteria(path)
original = bacteria(path)
globalNFE = 0
dAttr = 0.1
wAttr = 0.2
hRep = dAttr
wRep = 10

# Graficar
fitness_history = []
interaction_history = []
nfe_history = []
for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))


for _ in range(iteraciones):
    for b in poblacion:
        b.tumboNado(tumbo)
        b.autoEvalua()
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE
    best = max(poblacion, key=lambda x: x.fitness)
    if veryBest is None or best.fitness > veryBest.fitness:
        clonaBest(veryBest, best)

    fitness_history.append(veryBest.fitness)
    interaction_history.append(veryBest.interaction)
    nfe_history.append(globalNFE)

    print("Interacción:", veryBest.interaction, "Fitness:", veryBest.fitness, "NFE:", globalNFE, "Población:", len(poblacion))
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)


veryBest.showGenome()


plt.figure(figsize=(12, 6))
plt.plot(fitness_history, label="Fitness", marker='o')
plt.plot(interaction_history, label="Interacción", marker='o')
plt.plot(nfe_history, label="NFE Global", marker='o')
plt.xlabel("Iteraciones")
plt.ylabel("Valores")
plt.title("Evolución de Fitness, Interacción y NFE Global")
plt.legend()
plt.show()
