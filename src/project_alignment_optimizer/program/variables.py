import pathlib
# Variables
# ---------------------
# TODO: Refactor y pasar a un archivo las variables
class Var():
    def __init__(self):

        # Define la puntuacion a tener en cuenta para el calculo del score
        self._match = 1
        self._mismatch = -1
        self._gapPenalty = -1

        # Define el formato en el que recibiremos el archivo con el alineamiento inicial, siendo:
        # 0 -> fasta
        # 1 -> ???
        self._fileFormat = 0

        # Define el numero minimo de secuencias que tiene que tener el alineamiento
        self._nMinSequences = 50

        # Define la posicion en el array de la Secuencia Query
        self._querySequence = 0

        # Define de donde obtendremos las secuencias homologas, siendo:
        # 0 -> Base de datos local del usuario
        # 1 -> BLAST
        self._dbHomologousSequences = 0

        # Define de donde obtendremos las secuencias homologas de BLAST, siendo:
        # 0 -> corriendo de forma local
        # 1 -> corriendo desde la API ENTREZ
        self._dbBlast = 0

        self._clustalWPath = str(pathlib.Path(__file__).parent.parent.parent.parent.resolve()) + "/clustalw2"

        self._nToRemove = 20

    # ---------------------
    # Setters
    # ---------------------

    def setMatch(self,m):
        self._match = m

    def setMismatch(self,mm):
        self._mismatch = mm

    def setGapPenalty(self,gap):
        self._gapPenalty = gap

    def setFileFormat(self,ff):
        self._fileFormat = ff

    def setNMinSequences(self,nMin):
        self._nMinSequences = nMin

    def setQuerySequence(self,qs):
        self._querySequence = qs

    def setDbHomologousSequences(self,db):
        self._dbHomologousSequences = db

    def setDbBlast(self,db):
        self._dbBlast = db

    def setClustalWPath(self,aPath):
        self._clustalWPath = aPath
    
    def setNToRemove(self,n):
        self._nToRemove = n


    # ---------------------
    # Getters
    # ---------------------

    def match(self):
        return self._match

    def mismatch(self):
        return self._mismatch

    def gapPenalty(self):
        return self._gapPenalty

    def fileFormat(self):
        return self._fileFormat

    def nMinSequences(self):
        return self._nMinSequences

    def querySequence(self):
        return self._querySequence

    def dbHomologousSequences(self):
        return self._dbHomologousSequences

    def clustalWPath(self):
        return self._clustalWPath

    def dbBlast(self):
        return self._dbBlast

    def nToRemove(self):
        return self._nToRemove