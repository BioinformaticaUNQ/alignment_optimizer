import program.variables as var


# ---------------------
# Funciones Principales
# ---------------------

def configureVariables():
  # Pido por consola todas las variables a utilizar por el sistema
  printAndLog("Solicito Variables")

def loadFile():
  # Obtengo el archivo con el alineamiento inicial
  printAndLog("Solicito Variables")

def calculateScore(alignment):
  finalScore = 0
  for u in range(0,len(alignment)):
    for v in range(u+1,len(alignment)):
      finalScore += calculateScorePar(alignment[u],alignment[v])
  return finalScore

def filterAlignment(alignment):
  # Filtro el alineamiento pasado por parametro
  newAlignment = alignment
  return newAlignment

def getOriginalSequences(alignment):
  # Obtengo del alineamiento pasado por parametro las secuencias originales
  originalSequences = alignment
  return originalSequences

def generateAlignment(originalSequences):
  # Genero el nuevo alineamiento
  newAlignment = originalSequences
  return newAlignment

def generateTree(alignment):
  # Genero el árbol filogenético
  tree = alignment
  return tree


# ---------------------
# Funciones Auxiliates
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
    score += nucleotideCost(u[i],v[i])
  return score

def nucleotideCost(u,v):
    match = var.match()
    mismatch = var.mismatch()
    gap = var.gap()
    if (u == v):
        if(u == "_"):
            return 0
        return match
    else:
        if(u == "_" or v == "_"):
            return gap
        return mismatch

def printAndLog(msg):
    print(msg)
    # Guardo en el log el msg