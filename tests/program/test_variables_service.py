import pytest
import dotenv

from project_alignment_optimizer.program.constants import DB_BLAST, DB_HOMOLOGOUS_SEQUENCES, FILE_FORMAT, GAP_PENALTY, MATCH, MIN_SEQUENCES, MISMATCH, PURIFY_AMINO, RESET_VALUES
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
    assert variables_service.getVariableIntEnv(DB_BLAST) == dictVariablesValues[DB_BLAST]
    assert variables_service.getVariableIntEnv(PURIFY_AMINO) == dictVariablesValues[PURIFY_AMINO]

def test_getDictVariablesWithAllInfo():
    dictVariablesWithAllInfo = variables_service.getDictVariablesWithAllInfo()

    assert variables_service.getVariableIntEnv(MATCH) == dictVariablesWithAllInfo[MATCH]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(MISMATCH) == dictVariablesWithAllInfo[MISMATCH]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(GAP_PENALTY) == dictVariablesWithAllInfo[GAP_PENALTY]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(FILE_FORMAT) == dictVariablesWithAllInfo[FILE_FORMAT]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(MIN_SEQUENCES) == dictVariablesWithAllInfo[MIN_SEQUENCES]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES) == dictVariablesWithAllInfo[DB_HOMOLOGOUS_SEQUENCES]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(DB_BLAST) == dictVariablesWithAllInfo[DB_BLAST]['CURRENT_VALUE']
    assert variables_service.getVariableIntEnv(PURIFY_AMINO) == dictVariablesWithAllInfo[PURIFY_AMINO]['CURRENT_VALUE']

def test_tables_cant_be_tested():
    variables_service.getAllVariablesTable(SimpleNamespace(**{'file':None, 'query_sequence_header':None}))
    variables_service.getAllVariablesTableWithDescription()

def test_config_set_key_new_value_if():
    dotenv_file = dotenv.find_dotenv('config.env')
    new_match_value = -100

    old_match_value = variables_service.getVariableIntEnv(MATCH)
    assert new_match_value != old_match_value

    was_configured = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':MATCH, 'value':new_match_value}))
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables
    assert was_configured == True
    assert variables_service.getVariableIntEnv(MATCH) == new_match_value

    was_configured = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':MATCH, 'value':old_match_value}))
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables
    assert was_configured == True
    assert variables_service.getVariableIntEnv(MATCH) == old_match_value

def test_config_set_key_new_value_else():
    was_configured = variables_service.config_set_key_new_value(SimpleNamespace(**{'key':'ERROR', 'value':None}))
    assert was_configured == False

def test_resetDefaultValues():
    dotenv_file = dotenv.find_dotenv('config.env')

    old_match_value = variables_service.getVariableIntEnv(MATCH)
    old_mismatch_value = variables_service.getVariableIntEnv(MISMATCH)
    old_gap_penalty_value = variables_service.getVariableIntEnv(GAP_PENALTY)
    old_file_format_value = variables_service.getVariableIntEnv(FILE_FORMAT)
    old_min_sequences_value = variables_service.getVariableIntEnv(MIN_SEQUENCES)
    old_db_homologous_sequences_value = variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES)
    old_db_blast_value = variables_service.getVariableIntEnv(DB_BLAST)
    old_purify_amino_value = variables_service.getVariableIntEnv(PURIFY_AMINO)

    variables_service.resetDefaultValues()
    dotenv.load_dotenv(dotenv_file, override=True) # Need to reload de os.env variables

    assert variables_service.getVariableIntEnv(MATCH) == RESET_VALUES[MATCH]
    assert variables_service.getVariableIntEnv(MISMATCH) == RESET_VALUES[MISMATCH]
    assert variables_service.getVariableIntEnv(GAP_PENALTY) == RESET_VALUES[GAP_PENALTY]
    assert variables_service.getVariableIntEnv(FILE_FORMAT) == RESET_VALUES[FILE_FORMAT]
    assert variables_service.getVariableIntEnv(MIN_SEQUENCES) == RESET_VALUES[MIN_SEQUENCES]
    assert variables_service.getVariableIntEnv(DB_HOMOLOGOUS_SEQUENCES) == RESET_VALUES[DB_HOMOLOGOUS_SEQUENCES]
    assert variables_service.getVariableIntEnv(DB_BLAST) == RESET_VALUES[DB_BLAST]
    assert variables_service.getVariableIntEnv(PURIFY_AMINO) == RESET_VALUES[PURIFY_AMINO]

    # Set the variables 
    variables_service.setVariableEnv(MATCH, old_match_value)
    variables_service.setVariableEnv(MISMATCH, old_mismatch_value)
    variables_service.setVariableEnv(GAP_PENALTY, old_gap_penalty_value)
    variables_service.setVariableEnv(FILE_FORMAT, old_file_format_value)
    variables_service.setVariableEnv(MIN_SEQUENCES, old_min_sequences_value)
    variables_service.setVariableEnv(DB_HOMOLOGOUS_SEQUENCES, old_db_homologous_sequences_value)
    variables_service.setVariableEnv(DB_BLAST, old_db_blast_value)
    variables_service.setVariableEnv(PURIFY_AMINO, old_purify_amino_value)
