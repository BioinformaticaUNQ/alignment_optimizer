import os
import dotenv
from prettytable import PrettyTable
from project_alignment_optimizer.program.constants import ALL_ENV_VARIABLES

dotenv_file = dotenv.find_dotenv('config.env')
dotenv.load_dotenv(dotenv_file)

def getVariableIntEnv(key):
    return int(os.environ[key])

def getVariableStrEnv(key):
    return os.environ[key]

def getDictVariables(valuesAsInteger = False):
    dictVariables = {}
    # Recorre todas las variables del .env
    for variable_env in ALL_ENV_VARIABLES:
        if valuesAsInteger:
            # Agrego al dict con el value tipo int
            dictVariables[variable_env] = int(os.environ[variable_env])
        else:
            # Agrego al dict con el value tipo string
            dictVariables[variable_env] = os.environ[variable_env]
    return dictVariables

def setVariableEnv(key, value):
    dotenv.set_key(dotenv_file, key, str(value))

def getAllVariablesTable():
    table = PrettyTable(['Key', 'Value'])
    table.align = 'l' # Align a la izquierda l -> left
    dictVariables = getDictVariables()
    for variable_env in dictVariables.keys():
        table.add_row([variable_env, dictVariables[variable_env]])

    return table