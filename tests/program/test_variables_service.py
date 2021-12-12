import pytest
import dotenv

from project_alignment_optimizer.program.constants import ADMIT_HOMOLOGOUS, ADMIT_HOMOLOGOUS_TYPE, DB_HOMOLOGOUS_SEQUENCES, DB_HOMOLOGOUS_SEQUENCES_TYPE, FILE_FORMAT, FILE_FORMAT_TYPE, GAP_PENALTY, GAPEXT, GAPOPEN, MATCH, MATRIX, MIN_SEQUENCES, MISMATCH, N_HOMOLOGOUS_SEQUENCES, PURIFY_END, PURIFY_START, RESET_VALUES
from project_alignment_optimizer.program import variables_service
from types import SimpleNamespace


def test_getVariableIntEnv_and_setVariableEnv():
    dotenv_file = dotenv.find_dotenv('config.env')
    new_match_value = -100

    old_match_value = variables_service.getVariableIntEnv(MATCH)
    assert new_match_value != old_match_value

    variables_service.setVariableEnv(MATCH, new_match_value)
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables
    assert variables_service.getVariableIntEnv(MATCH) == new_match_value

    variables_service.setVariableEnv(MATCH, old_match_value)
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables
    assert variables_service.getVariableIntEnv(MATCH) == old_match_value

def test_getDictVariablesValues():
    dictVariablesValues = variables_service.getDictVariablesValues()

    assert variables_service.getVariableIntEnv(MATCH) == dictVariablesValues[MATCH]
    assert variables_service.getVariableIntEnv(MISMATCH) == dictVariablesValues[MISMATCH]
    assert variables_service.getVariableIntEnv(GAP_PENALTY) == dictVariablesValues[GAP_PENALTY]
    assert variables_service.getVariableIntEnv(FILE_FORMAT) == dictVariablesValues[FILE_FORMAT]
    assert variables_service.getVariableIntEnv(MIN_SEQUENCES) == dictVariablesValues[MIN_SEQUENCES]
    assert variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES) == dictVariablesValues[DB_HOMOLOGOUS_SEQUENCES]
    assert variables_service.getVariableIntEnv(PURIFY_START) == dictVariablesValues[PURIFY_START]
    assert variables_service.getVariableIntEnv(PURIFY_END) == dictVariablesValues[PURIFY_END]

def test_getDictVariablesWithAllInfo():
    dictVariablesWithAllInfo = variables_service.getDictVariablesWithAllInfo()

    assert variables_service.getVariableIntEnv(MATCH) == dictVariablesWithAllInfo[MATCH]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(MISMATCH) == dictVariablesWithAllInfo[MISMATCH]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(GAP_PENALTY) == dictVariablesWithAllInfo[GAP_PENALTY]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(FILE_FORMAT) == dictVariablesWithAllInfo[FILE_FORMAT]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(MIN_SEQUENCES) == dictVariablesWithAllInfo[MIN_SEQUENCES]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES) == dictVariablesWithAllInfo[DB_HOMOLOGOUS_SEQUENCES]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(PURIFY_START) == dictVariablesWithAllInfo[PURIFY_START]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(PURIFY_END) == dictVariablesWithAllInfo[PURIFY_END]['CURRENT_VALUE']

def test_tables_cant_be_tested():
    variables_service.getAllVariablesTable(SimpleNamespace(**{'file':None, 'query_sequence_header':None, 'homologous_sequences_path':None}))
    variables_service.getAllVariablesTable(SimpleNamespace(**{'file':None, 'query_sequence_header':None, 'homologous_sequences_path':'TEST'}))
    variables_service.getAllVariablesTableWithDescription()

def test_getDictTypes_OK():
    dictTypes = variables_service.getDictTypes(ADMIT_HOMOLOGOUS)
    assert ADMIT_HOMOLOGOUS_TYPE == dictTypes

    dictTypes = variables_service.getDictTypes(FILE_FORMAT)
    assert FILE_FORMAT_TYPE == dictTypes

    dictTypes = variables_service.getDictTypes(DB_HOMOLOGOUS_SEQUENCES)
    assert DB_HOMOLOGOUS_SEQUENCES_TYPE == dictTypes

def test_getDictTypes_None():
    key = 'NOT_TYPE_VALUE'
    dictTypes = variables_service.getDictTypes(key)
    assert dictTypes == None

def test_config_set_key_new_value_if():
    dotenv_file = dotenv.find_dotenv('config.env')
    new_match_value = '-100'

    old_match_value = str(variables_service.getVariableIntEnv(MATCH))
    assert int(new_match_value) != old_match_value

    was_configured, message = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':MATCH, 'value':new_match_value}))
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables
    assert was_configured == True
    assert variables_service.getVariableIntEnv(MATCH) == int(new_match_value)

    was_configured, message = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':MATCH, 'value':old_match_value}))
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables
    assert was_configured == True
    assert variables_service.getVariableIntEnv(MATCH) == int(old_match_value)

def test_config_set_key_new_value_else():
    was_configured, message = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':'ERROR', 'value':None}))
    assert was_configured == False

def test_config_set_key_new_value_else_else():
    was_configured, message = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':MATCH, 'value':'ERROR'}))
    assert was_configured == False

def test_resetDefaultValues():
    dotenv_file = dotenv.find_dotenv('config.env')

    old_match_value = variables_service.getVariableIntEnv(MATCH)
    old_mismatch_value = variables_service.getVariableIntEnv(MISMATCH)
    old_gap_penalty_value = variables_service.getVariableIntEnv(GAP_PENALTY)
    old_file_format_value = variables_service.getVariableIntEnv(FILE_FORMAT)
    old_min_sequences_value = variables_service.getVariableIntEnv(MIN_SEQUENCES)
    old_db_homologous_sequences_value = variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES)
    old_purify_start_value = variables_service.getVariableIntEnv(PURIFY_START)
    old_purify_end_value = variables_service.getVariableIntEnv(PURIFY_END)

    variables_service.resetDefaultValues()
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables

    assert variables_service.getVariableIntEnv(MATCH) == RESET_VALUES[MATCH]
    assert variables_service.getVariableIntEnv(MISMATCH) == RESET_VALUES[MISMATCH]
    assert variables_service.getVariableIntEnv(GAP_PENALTY) == RESET_VALUES[GAP_PENALTY]
    assert variables_service.getVariableIntEnv(FILE_FORMAT) == RESET_VALUES[FILE_FORMAT]
    assert variables_service.getVariableIntEnv(MIN_SEQUENCES) == RESET_VALUES[MIN_SEQUENCES]
    assert variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES) == RESET_VALUES[DB_HOMOLOGOUS_SEQUENCES]
    assert variables_service.getVariableIntEnv(PURIFY_START) == RESET_VALUES[PURIFY_START]
    assert variables_service.getVariableIntEnv(PURIFY_END) == RESET_VALUES[PURIFY_END]

    # Set the variables 
    variables_service.setVariableEnv(MATCH, old_match_value)
    variables_service.setVariableEnv(MISMATCH, old_mismatch_value)
    variables_service.setVariableEnv(GAP_PENALTY, old_gap_penalty_value)
    variables_service.setVariableEnv(FILE_FORMAT, old_file_format_value)
    variables_service.setVariableEnv(MIN_SEQUENCES, old_min_sequences_value)
    variables_service.setVariableEnv(DB_HOMOLOGOUS_SEQUENCES, old_db_homologous_sequences_value)
    variables_service.setVariableEnv(PURIFY_START, old_purify_start_value)
    variables_service.setVariableEnv(PURIFY_END, old_purify_end_value)

def test_validateKey_MATCH_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(MATCH, value)
    assert isValidInt == True

def test_validateKey_MATCH_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(MATCH, value)
    assert isValidInt == False

def test_validateKey_MISMATCH_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(MISMATCH, value)
    assert isValidInt == True

def test_validateKey_MISMATCH_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(MISMATCH, value)
    assert isValidInt == False

def test_validateKey_GAP_PENALTY_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(GAP_PENALTY, value)
    assert isValidInt == True

def test_validateKey_GAP_PENALTY_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(GAP_PENALTY, value)
    assert isValidInt == False

def test_validateKey_FILE_FORMAT_true():
    value = '0'
    isValidInt, message = variables_service.validateKey(FILE_FORMAT, value)
    assert isValidInt == True

def test_validateKey_FILE_FORMAT_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(FILE_FORMAT, value)
    assert isValidInt == False

    value = '999' # Out of range
    isValidInt, message = variables_service.validateKey(FILE_FORMAT, value)
    assert isValidInt == False

def test_validateKey_MIN_SEQUENCES_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(MIN_SEQUENCES, value)
    assert isValidInt == True

def test_validateKey_MIN_SEQUENCES_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(MIN_SEQUENCES, value)
    assert isValidInt == False

    value = '-10' # Negative integer
    isValidInt, message = variables_service.validateKey(MIN_SEQUENCES, value)
    assert isValidInt == False

def test_validateKey_DB_HOMOLOGOUS_SEQUENCES_true():
    value = '0'
    isValidInt, message = variables_service.validateKey(DB_HOMOLOGOUS_SEQUENCES, value)
    assert isValidInt == True

def test_validateKey_DB_HOMOLOGOUS_SEQUENCES_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(DB_HOMOLOGOUS_SEQUENCES, value)
    assert isValidInt == False

    value = '999' # Out of range
    isValidInt, message = variables_service.validateKey(DB_HOMOLOGOUS_SEQUENCES, value)
    assert isValidInt == False

def test_validateKey_N_HOMOLOGOUS_SEQUENCES_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(N_HOMOLOGOUS_SEQUENCES, value)
    assert isValidInt == True

def test_validateKey_N_HOMOLOGOUS_SEQUENCES_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(N_HOMOLOGOUS_SEQUENCES, value)
    assert isValidInt == False

    value = '-10' # Negative integer
    isValidInt, message = variables_service.validateKey(N_HOMOLOGOUS_SEQUENCES, value)
    assert isValidInt == False

def test_validateKey_ADMIT_HOMOLOGOUS_true():
    value = '0'
    isValidInt, message = variables_service.validateKey(ADMIT_HOMOLOGOUS, value)
    assert isValidInt == True

def test_validateKey_ADMIT_HOMOLOGOUS_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(ADMIT_HOMOLOGOUS, value)
    assert isValidInt == False

    value = '999' # Out of range
    isValidInt, message = variables_service.validateKey(ADMIT_HOMOLOGOUS, value)
    assert isValidInt == False

def test_validateKey_PURIFY_START_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(PURIFY_START, value)
    assert isValidInt == True

def test_validateKey_PURIFY_START_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(PURIFY_START, value)
    assert isValidInt == False

    value = '-10' # Negative integer
    isValidInt, message = variables_service.validateKey(PURIFY_START, value)
    assert isValidInt == False

def test_validateKey_PURIFY_END_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(PURIFY_END, value)
    assert isValidInt == True

def test_validateKey_PURIFY_END_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(PURIFY_END, value)
    assert isValidInt == False

    value = '-10' # Negative integer
    isValidInt, message = variables_service.validateKey(PURIFY_END, value)
    assert isValidInt == False

def test_validateKey_GAPOPEN_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(GAPOPEN, value)
    assert isValidInt == True

def test_validateKey_GAPOPEN_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(GAPOPEN, value)
    assert isValidInt == False

def test_validateKey_GAPEXT_true():
    value = '1'
    isValidInt, message = variables_service.validateKey(GAPEXT, value)
    assert isValidInt == True

def test_validateKey_GAPEXT_false():
    value = 'NOT_VALID_INT'
    isValidInt, message = variables_service.validateKey(GAPEXT, value)
    assert isValidInt == False

def test_validateKey_MATRIX():
    value = 'TEST'
    isValidStr, message = variables_service.validateKey(MATRIX, value)
    assert isValidStr == True

def test_validateKey_invalidKey():
    key = 'INVALID_KEY'
    value = 'TEST'
    isValidStr, message = variables_service.validateKey(key, value)
    assert isValidStr == False

def test_isValidHomologousSequencesPath_true():
    dotenv_file = dotenv.find_dotenv('config.env')
    
    old_admit_homologous_value = variables_service.getVariableIntEnv(ADMIT_HOMOLOGOUS)
    temp_admit_homologous_value = 1 # True

    variables_service.setVariableEnv(ADMIT_HOMOLOGOUS, temp_admit_homologous_value)
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables

    isValid = variables_service.isValidHomologousSequencesPath(SimpleNamespace(**{'file':None, 'query_sequence_header':None, 'homologous_sequences_path':'TEST'}))
    assert isValid == True

    variables_service.setVariableEnv(ADMIT_HOMOLOGOUS, old_admit_homologous_value)

def test_isValidHomologousSequencesPath_false():
    dotenv_file = dotenv.find_dotenv('config.env')
    
    old_admit_homologous_value = variables_service.getVariableIntEnv(ADMIT_HOMOLOGOUS)
    temp_not_admit_homologous_value = 0 # False

    variables_service.setVariableEnv(ADMIT_HOMOLOGOUS, temp_not_admit_homologous_value)
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables

    notIsValid = variables_service.isValidHomologousSequencesPath(SimpleNamespace(**{'file':None, 'query_sequence_header':None, 'homologous_sequences_path':'TEST'}))
    assert notIsValid == False

    variables_service.setVariableEnv(ADMIT_HOMOLOGOUS, old_admit_homologous_value)
