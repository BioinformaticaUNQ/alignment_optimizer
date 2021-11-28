# coding=utf-8
import pathlib
from variables import Var
import functions as func
import logging as log
import argparse
import sys


# ---------------------
# Logs
# ---------------------
log.basicConfig(filename='alignment_optimizer.log',
                format='%(asctime)s - %(levelname)s: %(message)s', level=log.DEBUG)
log.debug('Probando debug')
log.info('Probando info')
log.warning('Probando warning')
log.error('Probando error')


# ---------------------
# Programa Principal
# ---------------------

def main(args):

    # Inicializo las variables:
    var = Var()

    # Configuro las variables
    # --- Todo esto se lo tenemos que pedir al usuario y comentar cuales son los valores por defecto ---
    func.configureVariables()
    arg = getArgs(args)
    file = arg.file
    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    lastAlignment = func.loadFile(file)
    print(lastAlignment[3])

    if(len(lastAlignment)> var.nMinSequences()):
        # Obtengo las secuencias originales
        ungappedSequences = func.getungappedSequences(lastAlignment)
        print(lastAlignment[3])
        # Genero el alineamiento para obtener el score perteneciente al alineamiento inicial pasado por el usuario
        # Este alineamiento lo descarto, ya que no me sirve
        lastScore = func.generateAlignmentAndCalculateScore(ungappedSequences)

        # Realizo el filtrado de secuencias del alineamiento inicial pasado por el usuario
        aligmentFiltered = func.filterAlignment(lastAlignment)

        # Obtengo las secuencias originales
        ungappedSequences = func.getungappedSequences(aligmentFiltered)

        # Genero el nuevo alineamiento y calculo su score
        currentScore = func.generateAlignmentAndCalculateScore(ungappedSequences)
        currentAlignment = func.loadCurrentAlignment()

        # Genero nuevos alineamientos y sus scores correspondientes
        # mientras aumente el score actual o llegue al minimo de secuencias
        while(currentScore > lastScore and len(currentAlignment) > var.nMinSequences()):
            lastAlignment = currentAlignment
            lastScore = currentScore
            aligmentFiltered = func.filterAlignment(lastAlignment)
            ungappedSequences = func.getungappedSequences(aligmentFiltered)
            currentScore = func.generateAlignmentAndCalculateScore(ungappedSequences)
            currentAlignment = func.loadCurrentAlignment()

        print(len(currentAlignment))
        print(currentScore)
        # Genero el árbol filogenético y lo retorno
        tree = func.generateTree(currentAlignment)


        # TODO: Exportar el alineamiento final, a un path dado en los argumentos?
        # TODO: Hacer que el arbol se genere de verdad

        return tree
    else:
        func.printAndLog('Current amount of sequences provided is ' + str(len(lastAlignment)) + 
        ' and the minimum amount is ' + str(var.nMinSequences()))
       


def getArgs(args):
    parser = argparse.ArgumentParser(description='Alignment optimizer')
    # TODO: Agregar argumento de path de salida del arbol y del alineamiento final
    # TODO: Completar Help
    parser.add_argument('-f', '--file',
                           type=str,
                           help='File fasta format', required=False, default=(str(pathlib.Path().resolve()) + '/alignment.fasta'))
   # parser.add_argument('-f', '--file',
   #                     type=str,
   #                     help='File fasta format', required=False, default=(str(pathlib.Path().resolve()) + '/alintest1.fasta'))
    arg = parser.parse_args()
    return arg


if __name__ == '__main__':
    main(sys.argv[1:])

# Para correrlo:
# python src/project_alignment_optimizer/program/main.py
