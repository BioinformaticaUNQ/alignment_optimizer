from variables import Var
import functions as func
import logging as log

# ---------------------
# Logs
# ---------------------
log.basicConfig(filename='alignment_optimizer.log', format='%(asctime)s - %(levelname)s: %(message)s', level=log.DEBUG)
log.debug('Probando debug')
log.info('Probando info')
log.warning('Probando warning')
log.error('Probando error')

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

def main():

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
    nmin = Var().nMinSequences()
    while(currentScore > lastScore or len(currentAlignment) >= nmin):
        lastAlignment = currentAlignment
        lastScore = currentScore
        currentAlignment = func.filterAlignment(lastAlignment)
        originalSequences = func.getOriginalSequences(currentAlignment)
        currentAlignment = func.generateAlignment(originalSequences)
        currentScore = func.calculateScore(currentAlignment)

    # Genero el árbol filogenético y lo retorno
    tree = func.generateTree(currentAlignment)
    return tree


if __name__ == '__main__':
    main()

# Para correrlo:
# python src/project_alignment_optimizer/program/main.py