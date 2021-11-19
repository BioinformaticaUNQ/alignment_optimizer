# Variables
# ---------------------

# Define la puntuacion a tener en cuenta para el calculo del score
match = 1
mismatch = -1
gapPenalty = -1

# Define el formato en el que recibiremos el archivo con el alineamiento inicial, siendo:
# 0 -> fasta
# 1 -> ???
fileFormat = 0

# Define el numero minimo de secuencias que tiene que tener el alineamiento
nMinSequences = 25

# Define la posicion en el array de la Secuencia Query
querySequence = 0

# Define de donde obtendremos las secuencias homologas, siendo:
# 0 -> Base de datos local del usuario
# 1 -> BLAST
dbHomologousSequences = 0

# Define de donde obtendremos las secuencias homologas de BLAST, siendo:
# 0 -> corriendo de forma local
# 1 -> corriendo desde la API ENTREZ
dbBlast = 0


# ---------------------
# Setters
# ---------------------

def setMatch(m):
    match = m

def setMismatch(mm):
    mismatch = mm

def setGapPenalty(gap):
    gapPenalty = gap

def setFileFormat(ff):
    fileFormat = ff

def setNMinSequences(nMin):
    nMinSequences = nMin

def setQuerySequence(qs):
    querySequence = qs

def setDbHomologousSequences(db):
    dbHomologousSequences = db

def setDbBlast(db):
    dbBlast = db


# ---------------------
# Getters
# ---------------------

def match():
    return match

def mismatch():
    return mismatch

def gapPenalty():
    return gapPenalty

def fileFormat():
    return fileFormat

def nMinSequences():
    return nMinSequences

def querySequence():
    return querySequence

def dbHomologousSequences():
    return dbHomologousSequences

def dbBlast():
    return dbBlast