import os
import re
import dotenv
from prettytable import PrettyTable
from project_alignment_optimizer.program.constants import *

dotenv_file = dotenv.find_dotenv('config.env')
dotenv.load_dotenv(dotenv_file)

def getVariableIntEnv(key):
    return int(os.environ[key])

def getVariableFloatEnv(key):
    return float(os.environ[key])

def getVariableStrEnv(key):
    return os.environ[key]

def setVariableEnv(key, value):
    dotenv.set_key(dotenv_file, key, str(value))

def validateKey(key, value):
    if key == MATCH:
        return validateInt(value)
    elif key == MISMATCH:
        return validateInt(value)
    elif key == FILE_FORMAT:
        return validateTypeRange(key, value)
    elif key == MIN_SEQUENCES:
        return validatePositiveInt(value)
    elif key == DB_HOMOLOGOUS_SEQUENCES:
        return validateTypeRange(key, value)
    elif key == N_HOMOLOGOUS_SEQUENCES:
        return validatePositiveInt(value)
    elif key == ADMIT_HOMOLOGOUS:
        return validateTypeRange(key, value)
    elif key == PURIFY_START:
        return validatePositiveInt(value)
    elif key == PURIFY_END:
        return validatePositiveInt(value)
    elif key == GAPOPEN:
        return validateFloat(value)
    elif key == GAPEXT:
        return validateFloat(value)
    elif key == MATRIX:
        return True, None
    else:
        return False, 'Error'

def validateInt(value):
    if value == None or value == '':
        message = f"No new value was inserted."
        return False, message

    isInt = re.match("[-+]?\d+$", value)
    if not isInt:
        message = f"The new value: '{value}' is not valid as an integer."
        return False, message

    return True, None

def validateFloat(value):
    if value == None or value == '':
        message = f"No new value was inserted."
        return False, message

    isInt = re.match("[-+]?\d*\.?\d*$", value)
    if not isInt:
        message = f"The new value: '{value}' is not valid as a float."
        return False, message

    return True, None

def validatePositiveInt(value):
    if value == None or value == '':
        message = f"No new value was inserted."
        return False, message

    isInt, message = validateInt(value)
    if not isInt:
        return isInt, message

    ivalue = int(value)
    if ivalue < 0:
        message = f"The new value: '{ivalue}' is not a valid positive integer value."
        return False, message

    return True, None

def validateTypeRange(key, value):
    if value == None or value == '':
        message = f"No new value was inserted."
        return False, message

    isInt, message = validateInt(value)
    if not isInt:
        return isInt, message

    ivalue = int(value)
    dictTypes = getDictTypes(key)
    if not ivalue in dictTypes.keys():
        message = f"The new value: '{ivalue}' is not a valid type. ({dictTypes})"
        return False, message

    return True, None

def getDictTypes(key):
    for name, description, type in ALL_ENV_VARIABLES_WITH_DESCRIPTION:
        if name == key:
            return type

    return None

def getDictVariablesValues():
    dictVariablesValues = {}
    # Recorre todas las variables del .env
    for variable_env in ALL_ENV_VARIABLES_INT:
        # Agrego al dict con el value tipo int
        dictVariablesValues[variable_env] = int(os.environ[variable_env])

    for variable_env in ALL_ENV_VARIABLES_FLOATS:
        # Agrego al dict con el value tipo float
        dictVariablesValues[variable_env] = float(os.environ[variable_env])

    for variable_env in ALL_ENV_VARIABLES_STRINGS:
        # Agrego al dict con el value tipo str
        dictVariablesValues[variable_env] = os.environ[variable_env]

    return dictVariablesValues

def getDictVariablesWithAllInfo():
    dictVariablesInfo = {}
    # Recorre todas las variables del .env
    for name, description, type in ALL_ENV_VARIABLES_WITH_DESCRIPTION:
        tempDict = {}
        tempDict['CURRENT_VALUE'] = getValueVariable(name)
        tempDict['DEFAULT_VALUE'] = RESET_VALUES[name]
        tempDict['DESCRIPTION'] = description
        tempDict['TYPE'] = type

        dictVariablesInfo[name] = tempDict

    return dictVariablesInfo

def getValueVariable(variableName):
    if variableName in ALL_ENV_VARIABLES_STRINGS:
        return os.environ[variableName]
    if variableName in ALL_ENV_VARIABLES_INT:
        return int(os.environ[variableName])
    if variableName in ALL_ENV_VARIABLES_FLOATS:
        return float(os.environ[variableName])
    
    return None
    

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
        isValidKey, message = validateKey(key_upper, args.value)
        if isValidKey:
            setVariableEnv(key_upper, args.value)
            message = f"Key '{key_upper}' was modified correctly, new value = {args.value}."
            return True, message
        else:
            return False, message
    else:
        message= f"Key Unrecognized, valid keys:{ALL_ENV_VARIABLES}"
        return False, message

def isValidHomologousSequencesPath(args):
    admit_homologous = getVariableIntEnv(ADMIT_HOMOLOGOUS)
    
    if args.homologous_sequences_path is not None and admit_homologous != True:
        return False
    
    return True


def resetDefaultValues():
    setVariableEnv(MATCH, RESET_VALUES[MATCH])
    setVariableEnv(MISMATCH, RESET_VALUES[MISMATCH])
    setVariableEnv(FILE_FORMAT, RESET_VALUES[FILE_FORMAT])
    setVariableEnv(MIN_SEQUENCES, RESET_VALUES[MIN_SEQUENCES])
    setVariableEnv(DB_HOMOLOGOUS_SEQUENCES, RESET_VALUES[DB_HOMOLOGOUS_SEQUENCES])
    setVariableEnv(N_HOMOLOGOUS_SEQUENCES, RESET_VALUES[N_HOMOLOGOUS_SEQUENCES])
    setVariableEnv(ADMIT_HOMOLOGOUS, RESET_VALUES[ADMIT_HOMOLOGOUS])
    setVariableEnv(PURIFY_START, RESET_VALUES[PURIFY_START])
    setVariableEnv(PURIFY_END, RESET_VALUES[PURIFY_END])
    setVariableEnv(GAPOPEN, RESET_VALUES[GAPOPEN])
    setVariableEnv(GAPEXT, RESET_VALUES[GAPEXT])
    setVariableEnv(MATRIX, RESET_VALUES[MATRIX])
