import project_alignment_optimizer.program.variables as var
import logging as log


# ---------------------
# Funciones Principales
# ---------------------

def configureVariables():
  # Pido por consola todas las variables a utilizar por el sistema
  printAndLog("Solicito Variables")

def loadFile():
  # Obtengo el archivo con el alineamiento inicial
  # Validar el archivo
  printAndLog("Solicito el archivo")
  return []

def calculateScore(alignment):
  # Sacar el score del alineamiento obtenido por CLUSTAL
  finalScore = len(alignment)
  return finalScore

def filterAlignment(alignment):
  # Saco la secuencia que tiene mas gaps
  # Si hay mas de una, elimino la que tenga menor score con respecto a la secuencia query 
  # (uso alineamiento de a pares -> scorePar())
  # (Mas adelante, en esta funcion tambien usamos lo de las secuencias homologas)
  newAlignment = alignment
  return newAlignment

def getOriginalSequences(alignment):
  # Obtengo del alineamiento pasado por parametro las secuencias originales
  # Para eso elimino todos los gaps
  originalSequences = alignment
  return originalSequences

def generateAlignment(originalSequences):
  # Genero el nuevo alineamiento por medio de CLUSTAL
  newAlignment = originalSequences
  return newAlignment

def generateTree(alignment):
  # Genero el árbol filogenético
  tree = alignment
  return tree


# ---------------------
# Funciones Auxiliares
# ---------------------

def calculateScorePar(sequences):
  score = 0
  for u in range(0,len(sequences)):
    for v in range(u+1,len(sequences)):
      score += scorePar(sequences[u],sequences[v])
  return score

def scorePar(u,v):
  score = 0
  for i in range(0,len(u)):
    score += compareCost(u[i],v[i])
  return score

def compareCost(u,v):
    match = Var().match()
    mismatch = Var().mismatch()
    gap = Var().gapPenalty()
    if (u == v):
        if(u == "-"):
            return 0
        return match
    else:
        if(u == "-" or v == "-"):
            return gap
        return mismatch

def printAndLog(msg):
    print(msg)
    log.info(msg)