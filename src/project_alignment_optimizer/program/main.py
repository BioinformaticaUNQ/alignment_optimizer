# coding=utf-8
import pathlib
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

    # Configuro las variables
    # --- Todo esto se lo tenemos que pedir al usuario y comentar cuales son los valores por defecto ---
    func.configureVariables()
    arg = getArgs(args)
    file = arg.file
    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    currentAlignment = func.loadFile(file)

    if(func.alignmentHasNMinSequences(currentAlignment)):
        # Obtengo las secuencias originales
        ungappedSequences = func.getungappedSequences(currentAlignment)

        # Genero el alineamiento para obtener el score perteneciente al alineamiento inicial pasado por el usuario
        # Este alineamiento lo descarto, ya que no me sirve
        currentScore = func.generateAlignmentAndCalculateScore(ungappedSequences)

        # Genero nuevos alineamientos y sus scores correspondientes
        # Mientras aumente el score actual sigo aplicando los filtros
        bestAlignment = currentAlignment
        bestScore = currentScore
        better = True
        print("Current Alignment: " + str(len(currentAlignment)))
        while(better):
            lastAlignment = currentAlignment
            lastScore = currentScore
            # Hago el primer filtrado (Saco la secuencia que tenga mas aminoacidos en las columnas donde la query tenga gaps)
            print("FILTER 1")
            alignmentFiltered = func.filterAlignment(lastAlignment)
            currentAlignment , currentScore = generateNewAlignmentAndScore(alignmentFiltered)
            print("Current Alignment: " + str(len(currentAlignment)))
            if(currentScore > lastScore):
                bestAlignment = currentAlignment
                bestScore = currentScore
            else:
                # Hago el segundo filtrado (Saco la secuencia que tenga mas gaps de todo el alineamiento)
                print("FILTER 2")
                alignmentFiltered = func.filterAlignmentAlternative(lastAlignment)
                currentAlignment , currentScore = generateNewAlignmentAndScore(alignmentFiltered)
                print("Current Alignment: " + str(len(currentAlignment)))
                if(currentScore > lastScore):
                    bestAlignment = currentAlignment
                    bestScore = currentScore
                    print("MEJORO")
                else:
                    # TODO: Aqui se podria agregar un tercer filtrado (ver el agregado de homologas en otra situacion que no sea al llegar al nMin)
                    # Como no mejoro mas con ninguno de los filtrados termino con la busqueda
                    better = False
                    print("NO MEJORO")

        print(len(bestAlignment))
        print(bestScore)

        # Genero el árbol filogenético y lo retorno
        tree = func.generateTree(bestAlignment)

        # TODO: Exportar el alineamiento final, a un path dado en los argumentos?
        # TODO: Hacer que el arbol se genere de verdad

        return tree
    else:
        func.printAndLog('Current amount of sequences provided is ' + str(len(currentAlignment)) + 
        ' and it is less than the minimum of sequences established for the alignment')
       

def generateNewAlignmentAndScore(alignmentFiltered):
    ungappedSequences = func.getungappedSequences(alignmentFiltered)
    currentScore = func.generateAlignmentAndCalculateScore(ungappedSequences)
    currentAlignment = func.loadCurrentAlignment()
    return currentAlignment , currentScore


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
