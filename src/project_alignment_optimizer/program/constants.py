import pathlib

# Constantes ENV string
PATH = 'PATH'

MATCH = 'MATCH'
MATCH_DESCRIPTION = "Puntuacion entre dos caracteres iguales del alineamiento."

MISMATCH = 'MISMATCH'
MISMATCH_DESCRIPTION = "Puntuacion entre dos caracteres desiguales del alineamiento."

GAP_PENALTY = 'GAP_PENALTY'
GAP_PENALTY_DESCRIPTION = "Penalidad de gap."

FILE_FORMAT  = 'FILE_FORMAT'
FILE_FORMAT_DESCRIPTION = "Define el formato en el que recibiremos el archivo con el alineamiento inicial, siendo: 0 -> fasta, 1 -> ???."
FILE_FORMAT_TYPE = { 0:'fasta' }

MIN_SEQUENCES = 'MIN_SEQUENCES'
MIN_SEQUENCES_DESCRIPTION = "Define el numero minimo de secuencias que tiene que tener el alineamiento."

QUERY_SEQUENCE = 'QUERY_SEQUENCE'
QUERY_SEQUENCE_DESCRIPTION = "Define la posicion en el array de la Secuencia Query."

DB_HOMOLOGOUS_SEQUENCES = 'DB_HOMOLOGOUS_SEQUENCES'
DB_HOMOLOGOUS_SEQUENCES_DESCRIPTION = "Define de donde obtendremos las secuencias homologas, siendo: 0 -> base de datos local, 1 -> BLAST."
DB_HOMOLOGOUS_SEQUENCES_TYPE = { 0:'db local', 1:'BLAST'}

DB_BLAST = 'DB_BLAST'
DB_BLAST_DESCRIPTION = "Define de donde obtendremos las secuencias homologas de BLAST, siendo: 0 -> local, 1 -> API ENTREZ."
DB_BLAST_TYPE = { 0:'db local', 1:'API ENTREZ'}

PURIFY_AMINO = 'PURIFY_AMINO'
PURIFY_AMINO_DESCRIPTION = 'Numero de secuencias a eliminar del inicio y el final de cada alineamiento.'

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