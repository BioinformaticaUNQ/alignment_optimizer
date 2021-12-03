# coding=utf-8

from project_alignment_optimizer.program.constants import CLUSTALW_PATH
import logging as log
import pathlib
from Bio import SeqIO
from Bio import AlignIO
from Bio import Align
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
import re
import os
import sys

from Bio import Entrez,pairwise2
from Bio.SeqRecord import SeqRecord

import copy as c


# -------------------
# Funciones Principales
# ---------------------


def configureVariables():
    # Pido por consola todas las variables a utilizar por el sistema
    printAndLog("Solicito Variables")


def loadFile(filename):
    # Validar path
    # Validar formato
    # Validar que sea la cantidad de secuencias sea mayor al valor minimo
    # Obtengo el archivo con el alineamiento inicial
    checkIsValidPath(filename)
    with open(filename, "r") as handle:
        originalAlignment = AlignIO.read(handle, "fasta")
        if (any(originalAlignment)):
            printAndLog("Getting alignments")
            printAndLog(str(len(originalAlignment)) +
                        " aligned sequences were obtained correctly")
            return originalAlignment
        else:
            print("Invalid file format")


def checkIsValidPath(filename):
    if not os.path.isfile(filename):
        print("Invalid file path")
        sys.exit()


def calculateScore(alignment):
    # Sacar el score del alineamiento obtenido por CLUSTAL
    finalScore = len(alignment)
    return finalScore


def filterAlignment(anAlignment):
    lenSeq = len(anAlignment[0].seq)
    # Tomo el valor qe voy a cortar del inicio y del final para sacar los purificadores
    
    nToRemove = 20 # TODO: hacer este metodo  var.nToRemove
    # Agarro la sequencia query
    querySeq = anAlignment[0] # TODO: hacer este metodo   Var().querySequence()
    # Primero tengo que revisar que el largo de la cadena es mayor a las secuencias que voy a cortar del incio y el final
    if lenSeq > (nToRemove*2):
        start = nToRemove
        end = len(querySeq.seq)-nToRemove
    else:
        start = 0
        end = len(querySeq.seq)
    # Corto todas las secuencias del alineamiento por el numero obtenido en el paso anterior
    # Me fijo en que posiciones hay gaps y guardo los indices en una lista
    indices = []
    for i in range(start, end):
        if querySeq.seq[i] == '-':
            indices.append(i)
    # Recorro todas las secuencias y me voy quedando con los indices que esten en la lista
    seqToRemove = anAlignment[1]
    countSeqToRemove = -1
    for x in range(1,len(anAlignment)):
        seq = anAlignment[x].seq
        count = 0
        for y in indices:
            # Una vez obtenido esas secuencias filtradas, busco la que tenga mayor numero de aminoaciodos y esa es la que voy a eliminar
            if seq[y] != '-':
                count += 1
        if count > countSeqToRemove:
            seqToRemove = anAlignment[x]
    printAndLog("Sequence {id} has been removed from the alignment.".format(
        id=seqToRemove.id))
    sequences = list(filter(lambda seq: seq.id != seqToRemove.id, anAlignment))
    printAndLog(str(len(sequences)) + " sequences left.")
    print("Current Alignment: " + str(len(sequences)))
    response = checkToAdd(querySeq, sequences)
    print("Current Alignment: " + str(len(response)))
    return response


def filterAlignmentAlternative(anAlignment):
    # Saco la secuencia que tiene mas gaps
    querySeq = anAlignment[0] #TODO Var().querySequence() 
    mostGappedSeq = anAlignment[1]
    for ind in range(2,len(anAlignment)):
        mostGappedSeq = mostGapped(mostGappedSeq, anAlignment[ind], querySeq)
    printAndLog("Sequence {id} has been removed from the alignment.".format(
        id=mostGappedSeq.id))
    sequences = list(filter(lambda seq: seq.id != mostGappedSeq.id, anAlignment))
    printAndLog(str(len(sequences)) + " sequences left.")
    print("Current Alignment: " + str(len(sequences)))
    response = checkToAdd(querySeq, sequences)
    print("Current Alignment: " + str(len(response)))
    return response


def checkToAdd(seqQuery, sequences):
    # Por ahora agregamos para seguir manteniento el nMin
    # Mas adelante tendremos que agregar en otros momentos
    if len(sequences) < 50: #TODO: Var().nMinSequences()
        seqAux = selectSequenceToAdd(seqQuery, sequences)
        sequences.push(seqAux)
    return sequences


def selectSequenceToAdd(seqQuery, sequences):
    # Deberia devolver una secuencia que no este en el alineamiento y que no haya sido agregada previamente
    seqs = getHomologuesSequencesOrderedByMaxScore(seqQuery)
    # Duda: Trae todas las secuencias homologas que encuentra?
    # Duda: Podriamos hacerlo una sola vez, guardar las homologas e ir sacando las que ya agregamos?
    # Por ahora solo devuelvo la primera que devuelve, sin tener en cuenta lo de arriba
    for seq in seqs:
        if not seq in sequences:
            return seq


def alignmentHasNMinSequences(sequences):
    return len(sequences) >= 50 #TODO Var().nMinSequences()


def mostGapped(aSequence, anotherSequence, aQuerySequence):
    # Obtengo la secuencia con mayor cantidad de gaps
    if gaps(aSequence) > gaps(anotherSequence):
        return aSequence
    elif gaps(aSequence) == gaps(anotherSequence):
        aligner = Align.PairwiseAligner()
        alignment1 = aligner.align(aSequence.seq.ungap(), aQuerySequence.seq.ungap())
        alignment2 = aligner.align(anotherSequence.seq.ungap(),aQuerySequence.seq.ungap())
        if alignment1.score < alignment2.score:
            return aSequence
        else:
            return anotherSequence
    else:
        return anotherSequence


def gaps(aSequence):
    return len(aSequence.seq) - len(aSequence.seq.ungap())


def getungappedSequences(anAlignment):
    # Obtengo del alineamiento pasado por parametro las secuencias originales
    # Para eso elimino todos los gaps
    baseSequences = []
    ungappedSequences = c.deepcopy(anAlignment)
    for sequence in ungappedSequences:
        sequence.seq = sequence.seq.ungap()
        baseSequences.append(sequence)
    return baseSequences


def parseScore(aClustalOutputString):
    return re.search("(?<=Alignment Score )\d*", aClustalOutputString).group()


def generateAlignmentAndCalculateScore(originalSequences):
    # Genero el nuevo alineamiento por medio de CLUSTAL
    tempDir = str(pathlib.Path(__file__).parent.resolve())
    SeqIO.write(originalSequences, (tempDir + "/seqs.fasta"), "fasta")
    command = ClustalwCommandline(
        CLUSTALW_PATH, infile=(tempDir + "/seqs.fasta"))
    clusalAlignmentOutput = command()
    score = parseScore(clusalAlignmentOutput[0])
    printAndLog("New alignment Score: " + score)
    return score


def loadCurrentAlignment():
    tempDir = str(pathlib.Path(__file__).parent.resolve())
    return AlignIO.read(tempDir + "/seqs.aln", "clustal")


def generateTree(alignment):
    # Genero el árbol filogenético
    tree = alignment
    return tree


def getHomologuesSequencesOrderedByMaxScore(seqQuery):
    #Se pasa como parametro un SeqRecord
    #y devuelve una lista de SeqRecord 
    #ordenada de mayor a menor por score con la seqQuery
    #Habria que ver si solo traer, por ejemplo, un N definido de seq homologas
    result = []
    idsProteins = getIdsHomologuesSequences(seqQuery.id)
    if idsProteins:
        idString = ""
        for idProtein in idsProteins:
            idString = idString + idProtein + ","
        Entrez.email = "A.N.Other@example.com"
        idString = idString[0:len(idString)-1]
        handle = Entrez.efetch(db="protein", id=idString, rettype="gb", retmode="xml")
        output = Entrez.read(handle)
        result = getSequencesOrderedByMaxScore(seqQuery, output)
    return result


def calculateScorePar(sequences):
    score = 0
    for u in range(0, len(sequences)):
        for v in range(u+1, len(sequences)):
            score += scorePar(sequences[u], sequences[v])
    return score


def scorePar(u, v):
    score = 0
    for i in range(0, len(u)):
        score += compareCost(u[i], v[i])
    return score


def compareCost(u, v):
    match = 1  # var.match()
    mismatch = -1  # var.mismatch()
    gap = -1  # var.gapPenalty()
    if (u == v):
        return match
    else:
        if(u == "-" or v == "-"):
            return gap
        return mismatch


def printAndLog(msg):
    print(msg)
    log.info(msg)


def getIdsHomologuesSequences(idProtein):
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.efetch(db="protein", id=idProtein,rettype='ipg', retmode='xml')
    output = Entrez.read(handle)
    proteinList = output['IPGReport']['ProteinList']
    result = []
    if proteinList:
        for indx in range(0,len(proteinList)):
            result.append(proteinList[indx].attributes['accver'])
    return result


def getSequencesOrderedByMaxScore(seqQuery,output):
    result = []
    if output:
        for protein in output:
            if protein:
                descriptor = protein["GBSeq_locus"]
                seqHomologue = SeqRecord(Seq(protein["GBSeq_sequence"].upper()),id =descriptor,name=descriptor,description=descriptor)
                result.append((seqHomologue,pairwise2.align.globalxx(seqQuery.seq, seqHomologue.seq,score_only=True)))
        result.sort(key=takeSecond,reverse=True)
    return [i[0] for i in result] 


def takeSecond(elem):
    return elem[1]

