# coding=utf-8
import project_alignment_optimizer.program.variables_service as variables_service
import project_alignment_optimizer.program.functions as func
import logging as log
import argparse
import sys

from project_alignment_optimizer.program.constants import QUERY_RUN_ALIGN


# ---------------------
# Logs
# ---------------------
log.basicConfig(filename='alignment_optimizer.log', format='%(asctime)s - %(levelname)s: %(message)s', level=log.DEBUG)


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
                            help='Path to the align file',
                            required=True)

        parser.add_argument('-qs', '--query_sequence_header',
                            type=str,
                            help='Query sequence header name',
                            required=True)
        print(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        _align(args)

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

def _align(args):
    # Muestro la tabla de configuracion
    print(variables_service.getAllVariablesTable(args))
    
    # Pregunto si esta seguro que quiere ejecutar el comando con esa configuracion.
    if not func.query_yes_no(QUERY_RUN_ALIGN):
        return

    # Busco las variables globales
    env_variables = variables_service.getDictVariablesValues()

    func.align(args, env_variables)
   

def _config(args):
    if args.reset:
        variables_service.resetDefaultValues()
        print('File config.env restored to default values.')
    if args.key and args.value:
        was_configured = variables_service.config_set_key_new_value(args)
        if not was_configured:
            exit(1)


def _view_config():
    print(variables_service.getAllVariablesTableWithDescription())
       

if __name__ == '__main__':
    AlignmentOptimazer()
