# coding=utf-8
import pathlib
import project_alignment_optimizer.program.variables_service as variables_service
import functions as func
import logging as log
import argparse
import sys

from project_alignment_optimizer.program.constants import ALL_ENV_VARIABLES, MIN_SEQUENCES


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

    # Inicializo las variables:
    env_variables = variables_service.getDictVariables(True)

    # Configuro las variables
    # --- Todo esto se lo tenemos que pedir al usuario y comentar cuales son los valores por defecto ---
    func.configureVariables()

    file = args.file
    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    lastAlignment = func.loadFile(file)
    print(lastAlignment[3])

    if(len(lastAlignment)> env_variables[MIN_SEQUENCES]):
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
        while(currentScore > lastScore and len(currentAlignment) > env_variables[MIN_SEQUENCES]):
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
        ' and the minimum amount is ' + str(env_variables[MIN_SEQUENCES]))

def _config(args):
    if args.reset:
        pass # restart all variables.
    if args.key and args.value:
        key_upper = args.key.upper()
        dictVariables = variables_service.getDictVariables()
        if key_upper in dictVariables:
            variables_service.setVariableEnv(key_upper, args.value)
            print(f"Key '{key_upper}' was modified correctly, new value = {args.value}.")
        else:
            print(f'Key Unrecognized, valid keys:{ALL_ENV_VARIABLES}')
            exit(1)

def _view_config():
    print(variables_service.getAllVariablesTable())
       
if __name__ == '__main__':
    AlignmentOptimazer()

# Para correrlo:
# python src/project_alignment_optimizer/program/main.py