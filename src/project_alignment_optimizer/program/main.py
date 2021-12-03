# coding=utf-8
import pathlib
import project_alignment_optimizer.program.variables_service as variables_service
import functions as func
import logging as log
import argparse
import sys

from project_alignment_optimizer.program.constants import ALL_ENV_VARIABLES


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

class AlignmentOptimazer(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Description for AlignmentOptimazer',
            usage='''<command> [<args>]

            Commands:
            align           Run alignment code.
            config          Configuration management.
            view_config     View commands.
            ''')

        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def align(self):
        parser = argparse.ArgumentParser(description='Align .fasta file.')

        parser.add_argument('-f', '--file',
                            type=str,
                            help='File fasta format',
                            required=False,
                            default=(str(pathlib.Path().resolve()) + '/alignment.fasta'))

        args = parser.parse_args(sys.argv[2:])
        _aligh(args)

    def config(self):
        parser = argparse.ArgumentParser(description='Modify commands from the config.env file.')

        parser.add_argument('-r', '--reset',
                            type=bool,
                            help='Reset the config.env file to default.',
                            required=False)

        parser.add_argument('-k', '--key',
                            type=str,
                            help='Key from the config.env file.',
                            required=False)

        parser.add_argument('-v', '--value',
                            type=int,
                            help='Value for the key.',
                            required=False)

        args = parser.parse_args(sys.argv[2:])
        _config(args)
    
    def view_config(self):
        _view_config()


# ---------------------
# Funciones
# ---------------------

def _aligh(args):

    file = args.file
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
       
def _config(args):
    if args.reset:
        variables_service.resetDefaultValues()
        print('File config.env restored to default values.')
    if args.key and args.value:
        key_upper = args.key.upper()
        print(key_upper)
        dictVariables, dictDescriptions = variables_service.getDictVariables()
        if key_upper in dictVariables:
            variables_service.setVariableEnv(key_upper, args.value)
            print(f"Key '{key_upper}' was modified correctly, new value = {args.value}.")
        else:
            print(f'Key Unrecognized, valid keys:{ALL_ENV_VARIABLES}')
            exit(1)

def _view_config():
    print(variables_service.getAllVariablesTable())
       
def generateNewAlignmentAndScore(alignmentFiltered):
    ungappedSequences = func.getungappedSequences(alignmentFiltered)
    currentScore = func.generateAlignmentAndCalculateScore(ungappedSequences)
    currentAlignment = func.loadCurrentAlignment()
    return currentAlignment, currentScore

if __name__ == '__main__':
    AlignmentOptimazer()

# Para correrlo:
# python src/project_alignment_optimizer/program/main.py