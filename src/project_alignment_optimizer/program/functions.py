# coding=utf-8

from project_alignment_optimizer.program.constants import CLUSTALW_PATH, GAP_PENALTY, MATCH, MIN_SEQUENCES, MISMATCH, NOT_VALID_QUERY_NO, PURIFY_AMINO, VALID_QUERY_YES, DB_HOMOLOGOUS_SEQUENCES
from Bio import SeqIO, AlignIO, Align, Entrez, pairwise2
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
import logging as log
import pathlib
import re
import os
import sys
import copy as c

# ---------------------
# Funciones Principales
# ---------------------

def find_alignment_by_header(currentAlignment, query_sequence_header):
    for alignment in currentAlignment:
        if alignment.id == query_sequence_header:
            return alignment
    
    print(f"The header {query_sequence_header} has not been found.")
    sys.exit()

def align(args, env_variables):

    printAndLog("---------------------------------------")
    file = args.file
    query_sequence_header = args.query_sequence_header

    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    currentAlignment = loadFile(file)
    printAndLog("---------------------------------------")

    if(alignmentHasNMinSequences(currentAlignment, env_variables)):
        # Obtengo las secuencias homologas
        querySeq = find_alignment_by_header(currentAlignment, query_sequence_header)
        homologousSequences = getHomologousSequences(querySeq, currentAlignment, env_variables)
        printAndLog("---------------------------------------")

        # Obtengo las secuencias originales
        # TODO: Mejorar?
        seqsAux = c.deepcopy(currentAlignment)
        ungappedSequences = getungappedSequences(seqsAux)

        # Genero el alineamiento para obtener el score perteneciente al alineamiento inicial pasado por el usuario
        # Este alineamiento lo descarto, ya que no me sirve
        currentScore = generateAlignmentAndCalculateScore(ungappedSequences)

        # Genero nuevos alineamientos y sus scores correspondientes
        # Mientras aumente el score actual sigo aplicando los filtros
        bestAlignment = currentAlignment
        bestScore = currentScore
        better = True
        printAndLog("Current Alignment: " + str(len(currentAlignment)))
        while(better):
            lastAlignment = currentAlignment
            lastScore = currentScore
            printAndLog("---------------------------------------")
            # Hago el primer filtrado (Saco la secuencia que tenga mas aminoacidos en las columnas donde la query tenga gaps)
            printAndLog("FILTER 1 (Sequence that provides most gaps to the Query Sequence):")
            alignmentFiltered = filterSequenceThatProvidesMostGapsToQuery(lastAlignment, query_sequence_header, env_variables, homologousSequences)
            currentAlignment , currentScore = generateNewAlignmentAndScore(alignmentFiltered)
            printAndLog("Current Alignment: " + str(len(currentAlignment)))
            if(currentScore > lastScore):
                bestAlignment = currentAlignment
                bestScore = currentScore
                lastScore = currentScore
                lastAlignment = currentAlignment
                printAndLog("MEJORO 😁")
            else:
                printAndLog("NO MEJORO 😨")
                printAndLog("---------------------------------------")
                # Hago el segundo filtrado (Saco la secuencia que tenga mas gaps de todo el alineamiento)
                printAndLog("FILTER 2 (Sequence with most gaps):")
                alignmentFiltered = filterSequenceWithMostGaps(lastAlignment, query_sequence_header, env_variables, homologousSequences)
                currentAlignment , currentScore = generateNewAlignmentAndScore(alignmentFiltered)
                printAndLog("Current Alignment: " + str(len(currentAlignment)))
                if(currentScore > lastScore):
                    bestAlignment = currentAlignment
                    bestScore = currentScore
                    lastScore = currentScore
                    lastAlignment = currentAlignment
                    printAndLog("MEJORO 😀")
                else:
                    printAndLog("NO MEJORO 😨")
                    printAndLog("---------------------------------------")
                    # Hago el tercer filtrado (Agrego una secuencia homologa para ver si mejora el alineamiento)
                    printAndLog("FILTER 3 (Add a Homologous Sequence):")
                    if len(homologousSequences) > 0:
                        alignmentFiltered = addHomologousSequence(lastAlignment, homologousSequences)
                        currentAlignment , currentScore = generateNewAlignmentAndScore(alignmentFiltered)
                        printAndLog("Current Alignment: " + str(len(currentAlignment)))
                    else:
                        printAndLog("There are no more homologous sequences to add")
                        currentAlignment , currentScore = lastAlignment, lastScore
                    if(currentScore > lastScore):
                        bestAlignment = currentAlignment
                        bestScore = currentScore
                        lastScore = currentScore
                        lastAlignment = currentAlignment
                        printAndLog("MEJORO 😀")
                        printAndLog("---------------------------------------")
                    else:
                        # Como no mejoro mas con ninguno de los filtrados termino con la busqueda
                        better = False
                        printAndLog("NO MEJORO 😖")
                        printAndLog("---------------------------------------")

        printAndLog("The best alignment we obtained is " + str(len(bestAlignment)) + " sequences")
        printAndLog("The final score is " + str(bestScore))
        printAndLog("---------------------------------------")

        # Genero el árbol filogenético y lo retorno
        tree = generateTree(bestAlignment)
        printAndLog("Phylogenetic tree generated correctly")
        printAndLog("---------------------------------------")
        printAndLog("---------------------------------------")

        # TODO: Exportar el alineamiento final, a un path dado en los argumentos?
        # TODO: Hacer que el arbol se genere de verdad

        return tree
    else:
        printAndLog('Current amount of sequences provided is ' + str(len(currentAlignment)) + 
        ' and it is less than the minimum of sequences established for the alignment')


def loadFile(filename):
    # Validar path
    # Validar formato
    # Validar que sea la cantidad de secuencias sea mayor al valor minimo
    # Obtengo el archivo con el alineamiento inicial
    checkIsValidPath(filename)
    printAndLog("Getting Alignments")
    with open(filename, "r") as handle:
        originalAlignment = AlignIO.read(handle, "fasta")
        if (any(originalAlignment)):
            printAndLog(str(len(originalAlignment)) + " aligned sequences were obtained correctly")
            return originalAlignment
        else:
            print("Invalid file format")


def checkIsValidPath(filename):
    if not os.path.isfile(filename):
        print("Invalid file path")
        sys.exit()


def getHomologousSequences(querySeq, sequences, env_variables):
    printAndLog("Getting Homologous Sequences")
    db_hs = env_variables[DB_HOMOLOGOUS_SEQUENCES]
    if db_hs == 0:
        # TODO: Tendiamos que buscar en la base de datos local
        homologousSequences = getHomologousSequencesOrderedByMaxScore(querySeq)
    elif db_hs == 1:
        homologousSequences = getHomologousSequencesOrderedByMaxScore(querySeq)
    else:
        homologousSequences = getHomologousSequencesOrderedByMaxScore(querySeq)
    
    # TODO: Mejorar?
    response = list(filter(lambda seq: seq.id != sequences[0].id, homologousSequences))
    for ind in range(1,len(sequences)):
        response = list(filter(lambda seq: seq.id != sequences[ind].id, response))

    printAndLog(str(len(response)) + " homologous sequences were obtained correctly")
    return response


def filterSequenceThatProvidesMostGapsToQuery(anAlignment, query_sequence_header, env_variables, homologousSequences):
    min_sequences = env_variables[MIN_SEQUENCES]
    condition1 = len(anAlignment) > min_sequences
    condition2 = (len(anAlignment) <= min_sequences) and (len(homologousSequences) > 0)
    if condition1 or condition2:
        lenSeq = len(anAlignment[0].seq)
        # Tomo el valor que voy a cortar del inicio y del final para sacar los purificadores
        
        nToRemove = env_variables[PURIFY_AMINO]
        # Agarro la sequencia query
        querySeq = find_alignment_by_header(anAlignment, query_sequence_header)
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
                countSeqToRemove = count
                seqToRemove = anAlignment[x]
        printAndLog("Sequence {id} has been removed from the alignment.".format(id=seqToRemove.id))
        sequences = list(filter(lambda seq: seq.id != seqToRemove.id, anAlignment))
        printAndLog(str(len(sequences)) + " sequences left.")
        printAndLog("Current Alignment: " + str(len(sequences)))
        response = checkAndAddHomologousSequenceIfNeeds(sequences, homologousSequences, min_sequences)
        return response
    else:
        printAndLog("It cannot be filtered because there are no more homologous sequences to add and the minimum has already been reached")
        return anAlignment


def filterSequenceWithMostGaps(anAlignment, query_sequence_header, env_variables, homologousSequences):
    min_sequences = env_variables[MIN_SEQUENCES]
    condition1 = len(anAlignment) > min_sequences
    condition2 = (len(anAlignment) <= min_sequences) and (len(homologousSequences) > 0)
    if condition1 or condition2:
        # Saco la secuencia que tiene mas gaps
        querySeq = find_alignment_by_header(anAlignment, query_sequence_header)
        mostGappedSeq = anAlignment[1]

        for ind in range(2,len(anAlignment)):
            mostGappedSeq = mostGapped(mostGappedSeq, anAlignment[ind], querySeq)
        printAndLog("Sequence {id} has been removed from the alignment.".format(id=mostGappedSeq.id))
        sequences = list(filter(lambda seq: seq.id != mostGappedSeq.id, anAlignment))
        printAndLog(str(len(sequences)) + " sequences left.")
        printAndLog("Current Alignment: " + str(len(sequences)))
        response = checkAndAddHomologousSequenceIfNeeds(sequences, homologousSequences, min_sequences)
        return response
    else:
        printAndLog("It cannot be filtered because there are no more homologous sequences to add and the minimum has already been reached")
        return anAlignment


def checkAndAddHomologousSequenceIfNeeds(sequences, homologousSequences, min_sequences):
    if len(sequences) < min_sequences:
        response = addHomologousSequence(sequences, homologousSequences)
        return response
    else:
        return sequences


def addHomologousSequence(sequences, homologousSequences):
    # Agrega una secuencia que no este y que no haya sido agregada previamente
    seqToAdd = homologousSequences.pop(0)
    printAndLog(seqToAdd.id + " add to the Alignment")
    nNewSeq = len(seqToAdd.seq)
    nSeqAlignment = len(sequences[0].seq)
    nTotal = nSeqAlignment - nNewSeq
    seqFinal = seqToAdd + Seq("-"*nTotal)
    response = list(filter(lambda seq: True, sequences)) 
    # TODO: Mejorar?
    #response = sequences
    response.append(seqFinal)
    printAndLog("Current Alignment: " + str(len(response)))
    return response


def alignmentHasNMinSequences(sequences, env_variables):
    return len(sequences) >= env_variables[MIN_SEQUENCES]


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
    # TODO: Mejorar?
    baseSequences = []
    for sequence in anAlignment:
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


def getHomologousSequencesOrderedByMaxScore(seqQuery):
    #Se pasa como parametro un SeqRecord
    #y devuelve una lista de SeqRecord 
    #ordenada de mayor a menor por score con la seqQuery
    result = []
    idsProteins = getIdsHomologousSequences(seqQuery.id)
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


def printAndLog(msg):
    print(msg)
    log.info(msg)


def getIdsHomologousSequences(idProtein):
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
                seqHomologous = SeqRecord(Seq(protein["GBSeq_sequence"].upper()),id =descriptor,name=descriptor,description=descriptor)
                result.append((seqHomologous,pairwise2.align.globalxx(seqQuery.seq, seqHomologous.seq,score_only=True)))
        result.sort(key=takeSecond,reverse=True)
    return [i[0] for i in result] 


def takeSecond(elem):
    return elem[1]

def generateNewAlignmentAndScore(alignmentFiltered):
    ungappedSequences = getungappedSequences(alignmentFiltered)
    currentScore = generateAlignmentAndCalculateScore(ungappedSequences)
    currentAlignment = loadCurrentAlignment()
    return currentAlignment, currentScore

def query_yes_no(question):

    print(f"{question} (y/n):")
    input_value = input()

    if input_value is None:
        raise ValueError(f"invalid default answer: '{input_value}'")
    if input_value.lower() in VALID_QUERY_YES:
        return True
    elif input_value.lower() in NOT_VALID_QUERY_NO:
        return False
    else:
        print("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")
        query_yes_no(question)