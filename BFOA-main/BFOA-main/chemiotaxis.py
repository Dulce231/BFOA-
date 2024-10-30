import math
import random
import numpy
from bacteria import bacteria

class chemiotaxis():
    def __init__(self):
       self.parcialNFE = 0  
       self.d_attr = 0.1  
       self.w_attr = 0.2  
       self.h_rep = 0.1
       self.w_rep = 10 

    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 0.0
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore) ** 2.0
            total += d * math.exp(w * diff)
        return total

    def attract_repel(self, bacteria, poblacion):
        attract = self.compute_cell_interaction(bacteria, poblacion, -self.d_attr, -self.w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, self.h_rep, -self.w_rep)
        return attract + repel  

    def chemio(self, bacteria, poblacion):
        bacteria.interaction = self.attract_repel(bacteria, poblacion)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction

    def doChemioTaxis(self, poblacion):
        self.parcialNFE = 0
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion)
            self.parcialNFE += bacteria.NFE
            bacteria.NFE = 0

    def adaptarParametros(self, poblacion):
        # Calcula la desviación estándar del fitness de la población
        fitness_vals = [bacteria.fitness for bacteria in poblacion]
        std_fitness = numpy.std(fitness_vals)

        # Ajustar la dispersión del fitness
        if std_fitness < 5:
            self.d_attr *= 0.9
            self.h_rep *= 1.1
        else:
            self.d_attr *= 1.1
            self.h_rep *= 0.9
        # Ajustes límite para evitar valores extremos
        self.d_attr = min(max(self.d_attr, 0.05), 1.0)
        self.h_rep = min(max(self.h_rep, 0.05), 1.0)

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
        for i in range(num):
            poblacion.append(self.randomBacteria(path))
            poblacion.sort(key=lambda x: x.fitness)
            del poblacion[0]
