import os
import dotenv
from prettytable import PrettyTable
from project_alignment_optimizer.program.constants import ALL_ENV_VARIABLES, ALL_ENV_VARIABLES_WITH_DESCRIPTION, HOMOLOGOUS_SEQUENCES_PATH, MATCH, MISMATCH, GAP_PENALTY, FILE_FORMAT, MIN_SEQUENCES, DB_HOMOLOGOUS_SEQUENCES, PATH, QUERY_SEQUENCE_HEADER, RESET_VALUES, PATHLIB_ABSOLUTE

dotenv_file = dotenv.find_dotenv('config.env')
dotenv.load_dotenv(dotenv_file)

def getVariableIntEnv(key):
    return int(os.environ[key])

def setVariableEnv(key, value):
    # TODO FI: Validate value of the variable
    dotenv.set_key(dotenv_file, key, str(value))
    
def getDictVariablesValues():
    dictVariablesValues = {}
    # Recorre todas las variables del .env
    for variable_env in ALL_ENV_VARIABLES:
        # Agrego al dict con el value tipo int
        dictVariablesValues[variable_env] = int(os.environ[variable_env])
    return dictVariablesValues

def getDictVariablesWithAllInfo():
    dictVariablesInfo = {}
    # Recorre todas las variables del .env
    for name, description, type in ALL_ENV_VARIABLES_WITH_DESCRIPTION:
        tempDict = {}
        tempDict['CURRENT_VALUE'] = int(os.environ[name])
        tempDict['DEFAULT_VALUE'] = RESET_VALUES[name]
        tempDict['DESCRIPTION'] = description
        tempDict['TYPE'] = type

        dictVariablesInfo[name] = tempDict

    return dictVariablesInfo


def getAllVariablesTable(args):
    table = PrettyTable(['Key', 'Current value'])
    table.align = 'l' # Align a la izquierda l -> left
    dictVariables = getDictVariablesWithAllInfo()

    table.add_row([PATH, args.file])
    table.add_row([QUERY_SEQUENCE_HEADER, args.query_sequence_header])

    if args.homologous_sequences_path is not None:
        table.add_row([HOMOLOGOUS_SEQUENCES_PATH, PATHLIB_ABSOLUTE + args.homologous_sequences_path])
    else:
        table.add_row([HOMOLOGOUS_SEQUENCES_PATH, '-'])

    for variable_env_name in dictVariables.keys():
        if dictVariables[variable_env_name]['TYPE']:
            current_value = dictVariables[variable_env_name]['CURRENT_VALUE']
            type_info = dictVariables[variable_env_name]['TYPE'][current_value]
            current_value_with_type = f"{current_value} ({type_info})"
            table.add_row([ variable_env_name, current_value_with_type ])
        else:
            table.add_row([ variable_env_name, dictVariables[variable_env_name]['CURRENT_VALUE'] ])

    return table

def getAllVariablesTableWithDescription():
    table = PrettyTable(['Key', 'Current value', 'Default Value', 'Description'])
    table.align = 'l' # Align a la izquierda l -> left
    dictVariablesWithDescription = getDictVariablesWithAllInfo()
    for variable_env_name in dictVariablesWithDescription.keys():
        variable_info = dictVariablesWithDescription[variable_env_name]
        table.add_row([variable_env_name, variable_info['CURRENT_VALUE'],variable_info['DEFAULT_VALUE'] ,variable_info['DESCRIPTION']])

    return table

def config_set_key_new_value(args):
    key_upper = args.key.upper()
    if key_upper in ALL_ENV_VARIABLES:
        setVariableEnv(key_upper, args.value)
        print(f"Key '{key_upper}' was modified correctly, new value = {args.value}.")
        return True
    else:
        print(f'Key Unrecognized, valid keys:{ALL_ENV_VARIABLES}')
        return False

def resetDefaultValues():
    setVariableEnv(MATCH, RESET_VALUES[MATCH])
    setVariableEnv(MISMATCH, RESET_VALUES[MISMATCH])
    setVariableEnv(GAP_PENALTY, RESET_VALUES[GAP_PENALTY])
    setVariableEnv(FILE_FORMAT, RESET_VALUES[FILE_FORMAT])
    setVariableEnv(MIN_SEQUENCES, RESET_VALUES[MIN_SEQUENCES])
    setVariableEnv(DB_HOMOLOGOUS_SEQUENCES, RESET_VALUES[DB_HOMOLOGOUS_SEQUENCES])
    #setVariableEnv(PURIFY_AMINO, RESET_VALUES[PURIFY_AMINO])