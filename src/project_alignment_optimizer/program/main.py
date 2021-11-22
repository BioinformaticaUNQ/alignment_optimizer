from variables import Var
import functions as func
import logging as log
import argparse
import sys

# ---------------------
# Logs
# ---------------------
log.basicConfig(filename='alignment_optimizer.log', format='%(asctime)s - %(levelname)s: %(message)s', level=log.DEBUG)
log.debug('Probando debug')
log.info('Probando info')
log.warning('Probando warning')
log.error('Probando error')


# ---------------------
# Programa Principal
# ---------------------

def main(args):

    # Configuro las variables
    # --- Todo esto se lo tenemos que pedir al usuario y comentar cuales son los valores por defecto ---
    func.configureVariables()
    arg = getArgs(args)
    file = arg.file
    lastAlignment = func.loadFile(file)

    # Obtengo las secuencias originales
    originalSequences = func.getOriginalSequences(lastAlignment)

    # Genero el alineamiento para obtener el score perteneciente al alineamiento inicial pasado por el usuario
    # Este alineamiento lo descarto, ya que no me sirve
    lastAlignmentAux = func.generateAlignment(originalSequences)
    lastScore = func.calculateScore(lastAlignmentAux)

    # Realizo el filtrado de secuencias del alineamiento inicial pasado por el usuario
    aligmentFiltered = func.filterAlignment(lastAlignment)

    # Obtengo las secuencias originales
    originalSequences = func.getOriginalSequences(aligmentFiltered)

    # Genero el nuevo alineamiento y calculo su score
    currentAlignment = func.generateAlignment(originalSequences)
    currentScore = func.calculateScore(currentAlignment)

    # Genero nuevos alineamientos y sus scores correspondientes
    # mientras aumente el score actual o llegue al minimo de secuencias
    while(currentScore > lastScore or len(currentAlignment) >= Var().nMinSequences()):
        lastAlignment = currentAlignment
        lastScore = currentScore
        aligmentFiltered = func.filterAlignment(lastAlignment)
        originalSequences = func.getOriginalSequences(aligmentFiltered)
        currentAlignment = func.generateAlignment(originalSequences)
        currentScore = func.calculateScore(currentAlignment)

    # Genero el árbol filogenético y lo retorno
    tree = func.generateTree(currentAlignment)
    return tree

def getArgs(args):
    parser = argparse.ArgumentParser(description='Alignment optimizer')
    parser.add_argument('-f', '--file',
                    type=str,
                    help='File fasta format',required=True,default='')
    arg = parser.parse_args()
    return arg

if __name__ == '__main__':
   main(sys.argv[1:])

# Para correrlo:
# python src/project_alignment_optimizer/program/main.py