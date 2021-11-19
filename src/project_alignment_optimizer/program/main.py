import program.variables as var
import program.functions as func
import argparse


# ---------------------
# Variables
# ---------------------
lastScore = 0
currentScore = 0
lastAlignment = []
currenAlignment = []
originalSequences = []
tree = []


# ---------------------
# Programa Principal
# ---------------------

def main(args):
    args = argparse.parse_args(args)

    # Configuro las variables
    # --- Todo esto se lo tenemos que pedir al usuario y comentar cuales son los valores por defecto ---
    func.configureVariables()

    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    currentAlignment = func.loadFile()

    # Calculo el score del alineamiento inicial
    currentScore = func.calculateScore(currentAlignment)

    # Guardo el alineamiento inicial y su score como ultimo alineamiento
    lastAlignment = currentAlignment
    lastScore = currentScore

    # Realizo el filtrado de secuencias
    currentAlignment = func.filterAlignment(lastAlignment)

    # Obtengo las secuencias originales
    originalSequences = func.getOriginalSequences(currentAlignment)
   
    # Genero el nuevo alineamiento y calculo su score
    currentAlignment = func.generateAlignment(originalSequences)
    currentScore = func.calculateScore(currentAlignment)

    # Genero nuevos alineamientos y sus scores correspondientes 
    # mientras aumente el score actual o llegue al minimo de secuencias
    while(currentScore > lastScore or len(currentAlignment) >= var.nMinSequences()):
        lastAlignment = currentAlignment
        lastScore = currentScore
        currentAlignment = func.filterAlignment(lastAlignment)
        originalSequences = func.getOriginalSequences(currentAlignment)
        currentAlignment = func.generateAlignment(originalSequences)
        currentScore = func.calculateScore(currentAlignment)

    # Genero el árbol filogenético y lo retorno
    tree = func.generateTree(currentAlignment)
    return tree