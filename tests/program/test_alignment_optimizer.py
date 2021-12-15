import pytest
from project_alignment_optimizer.program import variables_service
from project_alignment_optimizer.program import functions
from project_alignment_optimizer.program.constants import MIN_SEQUENCES,PURIFY_START,PURIFY_END

def test_search_query_sequence():
    file = functions.loadFile('alignment.fasta')
    print('el header ' + '6QA2_A')
    query_seq = functions.find_alignment_by_header(file,'6QA2_A')
    print('la query_seq ' + str(query_seq))
    ungapped = query_seq.seq.ungap()
    assert ungapped == 'MTERTLVLIKPDGIERQLIGEIISRIERKGLTIAALQLRTVSAELASQHYAEHEGKPFFGSLLEFITSGPVVAAIVEGTAAIAAVRQLAGGTDPVQAAAPGTIRGDFALETQFNLVHGSDSAESAQREIALWFPGA'

def test_filter_aligment_with_more_amount_of_aacc():
   breakpoint()
   #args = Namespace(file='alignment.fasta', homologous_sequences_path=None, query_sequence_header='6QA2_A')
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   variables_service.setVariableEnv(MIN_SEQUENCES,90)
   mi_seq = variables_service.getVariableIntEnv(MIN_SEQUENCES)
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables,hom_path)
   #busco la que más gaps tiene
   sequence_with_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_seq)
   #valido que esa despues del primer filtrado no se encuentra 
   assert sequence_with_most_gaps.id == '1XQI_A'
   aligment_filtered = functions.filterSequence(0,currentAlignment, query_seq, homologousSequences, env_variables)
   filtered = list(filter(lambda al: al.seq == sequence_with_most_gaps.seq,aligment_filtered))
   notFiltered = list(filter(lambda al: al.seq == sequence_with_most_gaps.seq,currentAlignment))

   assert len(filtered) == 0 and len(notFiltered) == 1
  
def test_filter_sequence_provides_most_gapped_doesnt_inject_homologous():
   #no se agrega una homologa porque no llega al minimo de secuencias
   breakpoint()
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables,hom_path)
   seqs_homologous_before_load = homologousInCollection(currentAlignment,homologousSequences)

   assert len(seqs_homologous_before_load) == 0
   #busco la que más gaps tiene
   sequence_provides_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_seq)
   #aca ver que traiga las seqs no la que va a sacar 
   assert sequence_provides_most_gaps.id == '1XQI_A'

   sequences_without_seq_most_gaped = functions.filterSequence(0,currentAlignment, query_seq, homologousSequences, env_variables) 
   seqs_homologous_after_load = homologousInCollection(sequences_without_seq_most_gaped,homologousSequences)
   #ver que aca tendria que agregar la homologa .
   assert len(seqs_homologous_after_load)==0
   

def test_filter_sequence_most_gapped_inject_homologous():
   # se agrega una homologa porque llega al minimo de secuencias
   breakpoint()
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   variables_service.setVariableEnv(MIN_SEQUENCES,93)
   mi_seq = variables_service.getVariableIntEnv(MIN_SEQUENCES)
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables,hom_path)
   seqs_homologous_before_load= homologousInCollection(currentAlignment,homologousSequences)

   assert len(currentAlignment) == 93
   assert len(seqs_homologous_before_load)==0

   sequence_most_gap = functions.sequenceWithMostGaps(currentAlignment, query_seq,env_variables)
   assert sequence_most_gap.id == '4ANC_A'

   aligment_filtered =  functions.filterSequence(1,currentAlignment, query_seq, homologousSequences, env_variables)
   alignmentWithHomologousSequence, homologousSequence = functions.addHomologousSequence(
            aligment_filtered, homologousSequences)
   seqs_homologous_after_load= list(filter(lambda al: al.id == homologousSequence.id,alignmentWithHomologousSequence))
   #verificar que se agrego una secuencia homologa al alineamiento . y se mantiene la cantidad original de secuencias
   assert len(alignmentWithHomologousSequence) == 93
   assert len(seqs_homologous_after_load)==1

def test_trim_start_sequences():
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   breakpoint()
   variables_service.setVariableEnv(PURIFY_START,20)
   variables_service.setVariableEnv(PURIFY_END,0)
   p_s = variables_service.getVariableIntEnv(PURIFY_START)
   p_e = variables_service.getVariableIntEnv(PURIFY_END)
   env_variables = variables_service.getDictVariablesValues()

   sequence_no_trimmed = currentAlignment[0].seq
   assert len(sequence_no_trimmed) == 246
   breakpoint()
   aligment_trim = functions.trimPurifyingSequences(currentAlignment, env_variables)
   sequence_trimmed = aligment_trim[0].seq
   assert len(sequence_trimmed) == 226

def test_trim_start_and_end_sequences():
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)

   variables_service.setVariableEnv(PURIFY_START,20)
   variables_service.setVariableEnv(PURIFY_END,20)
   p_s = variables_service.getVariableIntEnv(PURIFY_START)
   p_e = variables_service.getVariableIntEnv(PURIFY_END)
   env_variables = variables_service.getDictVariablesValues()
   sequence_no_trimmed = currentAlignment[0].seq
   assert len(sequence_no_trimmed) == 246
   breakpoint()
   aligment_trim = functions.trimPurifyingSequences(currentAlignment, env_variables)
   sequence_trimmed = aligment_trim[0].seq
   assert len(sequence_trimmed) == 206


def test_default_trim_value_is_zero():
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)

   env_variables = variables_service.getDictVariablesValues()
   sequence_no_trimmed = currentAlignment[0].seq
   assert len(sequence_no_trimmed) == 246

   aligment_trim = functions.trimPurifyingSequences(currentAlignment, env_variables)
   sequence_trimmed = aligment_trim[0].seq
   assert len(sequence_no_trimmed) == 246

def test_after_first_filter_score_didn_improve_ramains_original_score():
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   variables_service.setVariableEnv(MIN_SEQUENCES,93)
   mi_seq = variables_service.getVariableIntEnv(MIN_SEQUENCES)

   ungappedSequences = functions.getUngappedSequences(currentAlignment)
   originalScore = functions.generateAlignmentAndCalculateScore(ungappedSequences, env_variables)

   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables,hom_path)
   #busco la que más gaps tiene
   breakpoint()
   sequence_with_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_seq)
   improve, alignmentFiltered, updatedScore =functions.executeFirstAlgorithm(currentAlignment, query_seq, homologousSequences, originalScore, env_variables)
   breakpoint()
   #valido que esa despues del primer filtrado no se encuentra 
   assert originalScore == updatedScore


def test_after_second_admit_homolougous_filter_score_improve_score():
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   hom_path = None
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   variables_service.setVariableEnv(MIN_SEQUENCES,93)
   mi_seq = variables_service.getVariableIntEnv(MIN_SEQUENCES)

   ungappedSequences = functions.getUngappedSequences(currentAlignment)
   originalScore = functions.generateAlignmentAndCalculateScore(ungappedSequences, env_variables)

   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables,hom_path)
   #busco la que más gaps tiene
   breakpoint()
   sequence_with_most_gaps = functions.sequenceWithMostGaps(currentAlignment, query_seq, env_variables)
   improve, alignmentFiltered, updatedScore =functions.executeSecondAlgorithm(currentAlignment, query_seq, homologousSequences, originalScore, env_variables)
   breakpoint()
   #valido que esa despues del primer filtrado no se encuentra 
   assert originalScore < updatedScore

def homologousInCollection(sequenceCollection,homologousSeqCollection):
   seq_hom_found =set()
   for sqh in homologousSeqCollection:
         seq_found = list(filter(lambda al: al.id == sqh.id,sequenceCollection))
         if len(seq_found)>0:
            seq_hom_found.add(set(seq_found))              
   return list(seq_hom_found)