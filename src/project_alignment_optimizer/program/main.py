# coding=utf-8
import project_alignment_optimizer.program.variables_service as variables_service
import project_alignment_optimizer.program.functions as func
import logging as log
import argparse
import sys

from project_alignment_optimizer.program.constants import QUERY_RUN_ALIGN, QUERY_RUN_RESET, TEMP_DIR


# ---------------------
# Logs
# ---------------------
log.basicConfig(filename=TEMP_DIR + '/alignment_optimizer.log', format='%(asctime)s - %(levelname)s: %(message)s', level=log.DEBUG)


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
 reset           Reset the config.env file to default.
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
                            help='Path to the align file',
                            required=True)

        parser.add_argument('-qsh', '--query_sequence_header',
                            type=str,
                            help='Query sequence header name',
                            required=True)
        
        parser.add_argument('-hsp', '--homologous_sequences_path',
                            type=str,
                            help='Path to the homologous sequences .fasta file (requires ADMIT_HOMOLOGOUS parameter enabled)',
                            required=False)

        args = parser.parse_args(sys.argv[2:])

        _align(args)

    def config(self):
        parser = argparse.ArgumentParser(description='Modify commands from the config.env file.')

        parser.add_argument('-k', '--key',
                            type=str,
                            help='Key from the config.env file.',
                            required=True)

        parser.add_argument('-v', '--value',
                            help='Value for the key.',
                            required=True)

        args = parser.parse_args(sys.argv[2:])
        _config(args)
    
    def view_config(self):
        _view_config()

    def reset(self):
        _reset_config()


# ---------------------
# Funciones
# ---------------------

def _align(args):
    # Valido que si se le pasa un path de --homologous_sequences_path debe de tener ENABLE el parametro ADMIT_HOMOLOGOUS
    if not variables_service.isValidHomologousSequencesPath(args):
        print("If the argument '--homologous_sequences_path' is passed, the 'ADMIT_HOMOLOGOUS' parameter must be enabled.")
        return

    # Muestro la tabla de configuracion
    print(variables_service.getAllVariablesTable(args))
    
    # Pregunto si esta seguro que quiere ejecutar el comando con esa configuracion.
    if not func.query_yes_no(QUERY_RUN_ALIGN):
        return

    # Busco las variables globales
    env_variables = variables_service.getDictVariablesValues()
    log.info("The program will run with the following parameters: \n" + str(variables_service.getAllVariablesTable(args)))

    func.align(args, env_variables)
   

def _config(args):
    was_configured, message = variables_service.config_set_key_new_value(args)
    print(message)
    if not was_configured:
        exit(1)


def _view_config():
    print(variables_service.getAllVariablesTableWithDescription())

def _reset_config():
    # Pregunto si esta seguro que quiere ejecutar el comando reset.
    if not func.query_yes_no(QUERY_RUN_RESET):
        return

    variables_service.resetDefaultValues()
    print('File config.env restored to default values.')
       

if __name__ == '__main__':
    AlignmentOptimazer()
