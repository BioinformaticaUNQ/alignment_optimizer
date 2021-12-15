import pathlib

# Constantes ENV string
PATH = 'PATH'

MATCH = 'MATCH'
MATCH_DESCRIPTION = "Defines the score given for each character match in the alignment."

MISMATCH = 'MISMATCH'
MISMATCH_DESCRIPTION = "Defines the score given for each character mismatch in the alignment."

QUERY_SEQUENCE_HEADER = 'QUERY_SEQUENCE_HEADER'
QUERY_SEQUENCE_HEADER_DESCRIPTION = "Defines the header to be used as query sequence when filtering."

FILE_FORMAT = 'FILE_FORMAT'
FILE_FORMAT_DESCRIPTION = "Selects the file format of the input file. Currently only supports 0 (FASTA)."
FILE_FORMAT_TYPE = {0: 'fasta'}

MIN_SEQUENCES = 'MIN_SEQUENCES'
MIN_SEQUENCES_DESCRIPTION = "Defines the minimun number of sequences inside the alignment."

DB_HOMOLOGOUS_SEQUENCES = 'DB_HOMOLOGOUS_SEQUENCES'
DB_HOMOLOGOUS_SEQUENCES_DESCRIPTION = "Selects where to look for homologous sequences, where: 0 -> ENTREZ API, 1 -> local database."
DB_HOMOLOGOUS_SEQUENCES_TYPE = {0: 'ENTREZ API', 1: 'local database'}

MATRIX = "MATRIX"
MATRIX_DESCRIPTION = "Selects Clustal's protein wheight matrix (BLOSUM, PAM, GONNET, ID or filename)."

GAPOPEN = "GAPOPEN"
GAPOPEN_DESCRIPTION = "Defines Clustal's gap opening penalty."

GAPEXT = "GAPEXT"
GAPEXT_DESCRIPTION = "Defines Clustal's gap extension penalty."

HOMOLOGOUS_SEQUENCES_PATH = 'HOMOLOGOUS_SEQUENCES_PATH'
HOMOLOGOUS_SEQUENCES_PATH_DESCRIPTION = "Defines the path where the homologous sequences are found."

N_HOMOLOGOUS_SEQUENCES = "N_HOMOLOGOUS_SEQUENCES"
N_HOMOLOGOUS_SEQUENCES_DESCRIPTION = "Defines the number of homologous sequences to use."

ADMIT_HOMOLOGOUS = "ADMIT_HOMOLOGOUS"
ADMIT_HOMOLOGOUS_DESCRIPTION = "Defines if agree to add homologous sequences to the alignment, where: 0 -> NOT Admit, 1 -> Admit."
ADMIT_HOMOLOGOUS_TYPE = {0: 'No', 1: 'Yes'}

PURIFY_START = 'PURIFY_START'
PURIFY_START_DESCRIPTION = "Defines how many aminoacids should be trimmed from the beginning of each sequence."

PURIFY_END = 'PURIFY_END'
PURIFY_END_DESCRIPTION = "Defines how many aminoacids should be trimmed at the end of each sequence."

ALL_ENV_VARIABLES = [ADMIT_HOMOLOGOUS, MATCH, MISMATCH, FILE_FORMAT, MIN_SEQUENCES,
                     DB_HOMOLOGOUS_SEQUENCES, N_HOMOLOGOUS_SEQUENCES, PURIFY_START, PURIFY_END, GAPOPEN, GAPEXT, MATRIX]
ALL_ENV_VARIABLES_INT = [ADMIT_HOMOLOGOUS, MATCH, MISMATCH, FILE_FORMAT, MIN_SEQUENCES,
                         DB_HOMOLOGOUS_SEQUENCES, N_HOMOLOGOUS_SEQUENCES, PURIFY_START, PURIFY_END]
ALL_ENV_VARIABLES_STRINGS = [MATRIX]
ALL_ENV_VARIABLES_FLOATS = [GAPOPEN, GAPEXT]

ALL_ENV_VARIABLES_WITH_DESCRIPTION = [
    (ADMIT_HOMOLOGOUS, ADMIT_HOMOLOGOUS_DESCRIPTION, ADMIT_HOMOLOGOUS_TYPE),
    (MATCH, MATCH_DESCRIPTION, None),
    (MISMATCH, MISMATCH_DESCRIPTION, None),
    (MIN_SEQUENCES, MIN_SEQUENCES_DESCRIPTION, None),
    (DB_HOMOLOGOUS_SEQUENCES, DB_HOMOLOGOUS_SEQUENCES_DESCRIPTION,
     DB_HOMOLOGOUS_SEQUENCES_TYPE),
    (FILE_FORMAT, FILE_FORMAT_DESCRIPTION,
     FILE_FORMAT_TYPE),
    (N_HOMOLOGOUS_SEQUENCES,
     N_HOMOLOGOUS_SEQUENCES_DESCRIPTION, None),
    (PURIFY_START, PURIFY_START_DESCRIPTION, None),
    (PURIFY_END, PURIFY_END_DESCRIPTION, None),
    (GAPOPEN, GAPOPEN_DESCRIPTION, None),
    (GAPEXT, GAPEXT_DESCRIPTION, None),
    (MATRIX, MATRIX_DESCRIPTION, None)
]

# Constants PATH Clustal
CLUSTALW_PATH = str(pathlib.Path(
    __file__).parent.parent.parent.parent.resolve()) + "/clustalw2"
PATHLIB_ABSOLUTE = str(pathlib.Path(__file__).parent.absolute())

# Constants temp files dir
TEMP_DIR = str(pathlib.Path(__file__).parent.resolve()) + "/temp"

# Constants query
QUERY_RUN_ALIGN = 'Are you sure you want to run the alignment with these settings?'
QUERY_RUN_RESET = 'Are you sure you want to run the reset command?'
VALID_QUERY_YES = {"yes": True, "y": True, '': True}
NOT_VALID_QUERY_NO = {"no": False, "n": False}

# Reset values
RESET_VALUES = {
    MATCH: 1,
    MISMATCH: -1,
    FILE_FORMAT: 0,
    MIN_SEQUENCES: 50,
    DB_HOMOLOGOUS_SEQUENCES: 0,
    N_HOMOLOGOUS_SEQUENCES: 20,
    ADMIT_HOMOLOGOUS: 0,
    PURIFY_START: 0,
    PURIFY_END: 0,
    GAPOPEN: 10.0,
    GAPEXT: 0.2,
    MATRIX: 'PAM'
}
