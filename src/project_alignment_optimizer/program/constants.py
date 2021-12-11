import pathlib

# Constantes ENV string
PATH = 'PATH'

MATCH = 'MATCH'
MATCH_DESCRIPTION = "Defines the score given for each character match in the alignment."

MISMATCH = 'MISMATCH'
MISMATCH_DESCRIPTION = "Defines the score given for each character mismatch in the alignment."

GAP_PENALTY = 'GAP_PENALTY'
GAP_PENALTY_DESCRIPTION = "Defines the penalty assigned for each gap opening in the alignment."
###################################################

QUERY_SEQUENCE_HEADER = 'QUERY_SEQUENCE_HEADER'
QUERY_SEQUENCE_HEADER_DESCRIPTION = "Defines the header to be used as query sequence when filtering."

FILE_FORMAT  = 'FILE_FORMAT'
FILE_FORMAT_DESCRIPTION = "Selects the file format of the input file. Currently only supports 0 (FASTA)."
FILE_FORMAT_TYPE = { 0:'fasta' }

MIN_SEQUENCES = 'MIN_SEQUENCES'
MIN_SEQUENCES_DESCRIPTION = "Defines the minimun number of sequences inside the alignment."

DB_HOMOLOGOUS_SEQUENCES = 'DB_HOMOLOGOUS_SEQUENCES'
DB_HOMOLOGOUS_SEQUENCES_DESCRIPTION = "Selects where to look for homologous sequences, where: 0 -> local database, 1 -> BLAST."
DB_HOMOLOGOUS_SEQUENCES_TYPE = { 0:'db local', 1:'BLAST'}
DB_HOMOLOGOUS_SEQUENCES_PATH = str(pathlib.Path(__file__).parent.absolute()) + "/dbLocal.fasta"
N_HOMOLOGOUS_SEQUENCES = 20

DB_BLAST = 'DB_BLAST'
DB_BLAST_DESCRIPTION = "Selects where to look for BLAST homologous sequences, where: 0 -> local BLAST, 1 -> ENTREZ API."
DB_BLAST_TYPE = { 0:'db local', 1:'API ENTREZ'}

PURIFY_AMINO = 'PURIFY_AMINO'
PURIFY_AMINO_DESCRIPTION = "Defines how many aminoacids should be trimmed from the beginning and end of each sequence."

ALL_ENV_VARIABLES = [MATCH, MISMATCH, GAP_PENALTY, FILE_FORMAT, MIN_SEQUENCES, DB_HOMOLOGOUS_SEQUENCES, DB_BLAST, PURIFY_AMINO]
ALL_ENV_VARIABLES_WITH_DESCRIPTION = [
                                        (MATCH, MATCH_DESCRIPTION, None),
                                        (MISMATCH, MISMATCH_DESCRIPTION, None),
                                        (GAP_PENALTY, GAP_PENALTY_DESCRIPTION, None),
                                        (MIN_SEQUENCES, MIN_SEQUENCES_DESCRIPTION, None),
                                        (PURIFY_AMINO, PURIFY_AMINO_DESCRIPTION, None),
                                        (DB_HOMOLOGOUS_SEQUENCES, DB_HOMOLOGOUS_SEQUENCES_DESCRIPTION, DB_HOMOLOGOUS_SEQUENCES_TYPE),
                                        (DB_BLAST, DB_BLAST_DESCRIPTION, DB_BLAST_TYPE),
                                        (FILE_FORMAT, FILE_FORMAT_DESCRIPTION, FILE_FORMAT_TYPE)
                                    ]

# Constants PATH Clustal
CLUSTALW_PATH = str(pathlib.Path(__file__).parent.parent.parent.parent.resolve()) + "/clustalw2"

# Constants query
QUERY_RUN_ALIGN = 'Are you sure you want to run the alignment with these settings?'
VALID_QUERY_YES = {"yes": True, "y": True, '': True}
NOT_VALID_QUERY_NO = {"no": False, "n": False}

# Reset values
RESET_VALUES = {
    MATCH: 1,
    MISMATCH: -1,
    GAP_PENALTY: -1,
    FILE_FORMAT: 0,
    MIN_SEQUENCES: 50,
    DB_HOMOLOGOUS_SEQUENCES: 0,
    DB_BLAST: 0,
    PURIFY_AMINO: 0
}