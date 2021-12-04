import os
import dotenv
from prettytable import PrettyTable
from project_alignment_optimizer.program.constants import ALL_ENV_VARIABLES_WITH_DESCRIPTION, MATCH, MISMATCH, GAP_PENALTY, FILE_FORMAT, MIN_SEQUENCES, PURIFY_AMINO, QUERY_SEQUENCE, DB_HOMOLOGOUS_SEQUENCES, DB_BLAST

dotenv_file = dotenv.find_dotenv('config.env')
dotenv.load_dotenv(dotenv_file)

def getVariableIntEnv(key):
    return int(os.environ[key])

def getVariableStrEnv(key):
    return os.environ[key]

def getDictVariables(valuesAsInteger = False):
    dictVariables = {}
    dictDescriptions = {}
    # Recorre todas las variables del .env
    for variable_env, variable_env_description in ALL_ENV_VARIABLES_WITH_DESCRIPTION:
        if valuesAsInteger:
            # Agrego al dict con el value tipo int
            dictVariables[variable_env] = int(os.environ[variable_env])
        else:
            # Agrego al dict con el value tipo string
            dictVariables[variable_env] = os.environ[variable_env]
        dictDescriptions[variable_env] = variable_env_description
    return dictVariables, dictDescriptions

def setVariableEnv(key, value):
    dotenv.set_key(dotenv_file, key, str(value))

def getAllVariablesTable():
    table = PrettyTable(['Key', 'Value', 'Description'])
    table.align = 'l' # Align a la izquierda l -> left
    dictVariables, dictDescriptions = getDictVariables()
    for variable_env in dictVariables.keys():
        table.add_row([variable_env, dictVariables[variable_env], dictDescriptions[variable_env]])

    return table

def resetDefaultValues():
    setVariableEnv(MATCH, 1)
    setVariableEnv(MISMATCH, -1)
    setVariableEnv(GAP_PENALTY, -1)
    setVariableEnv(FILE_FORMAT, -1)
    setVariableEnv(MIN_SEQUENCES, 35)
    setVariableEnv(QUERY_SEQUENCE, 0)
    setVariableEnv(DB_HOMOLOGOUS_SEQUENCES, 0)
    setVariableEnv(DB_BLAST, 0)
    setVariableEnv(PURIFY_AMINO, 20)