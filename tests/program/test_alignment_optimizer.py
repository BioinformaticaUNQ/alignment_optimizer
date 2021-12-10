import pytest
from project_alignment_optimizer.program import variables_service
from project_alignment_optimizer.program import functions
from project_alignment_optimizer.program.constants import MIN_SEQUENCES

def test_search_query_sequence():
    file = functions.loadFile('alignment.fasta')
    print('el header ' + '6QA2_A')
    query_seq = functions.find_alignment_by_header(file,'6QA2_A')
    print('la query_seq ' + str(query_seq))
    ungapped = query_seq.seq.ungap()
    assert ungapped == 'MTERTLVLIKPDGIERQLIGEIISRIERKGLTIAALQLRTVSAELASQHYAEHEGKPFFGSLLEFITSGPVVAAIVEGTAAIAAVRQLAGGTDPVQAAAPGTIRGDFALETQFNLVHGSDSAESAQREIALWFPGA'

def test_filter_aligment_with_more_amount_of_aacc():
   breakpoint()
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables)
   #busco la que más gaps tiene
   sequence_with_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_sec_header, env_variables)
   #valido que esa despues del primer filtrado no se encuentra 
   breakpoint()
   aligment_filtered = functions.filterSequenceThatProvidesMostGapsToQuery(currentAlignment, '6QA2_A', env_variables, homologousSequences)
   breakpoint()
   filtered = list(filter(lambda al: al.seq == sequence_with_most_gaps.seq,aligment_filtered))
   notFiltered = list(filter(lambda al: al.seq == sequence_with_most_gaps.seq,currentAlignment))

   assert len(filtered) == 0 and len(notFiltered == 1)
  
def test_filter_sequence_most_gapped_doesnt_inject_homologous():
   #no se agrega una homologa porque no llega al minimo de secuencias
   breakpoint()
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables)
   seqs_homologous_before_load= homologousInCollection(currentAlignment,homologousSequences)

   assert len(seqs_homologous_before_load)==0
   #busco la que más gaps tiene
   sequence_provides_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_sec_header, env_variables)
   #aca ver que traiga las seqs no la que va a sacar 
   breakpoint()
   aligment_filtered = functions.filterSequenceThatProvidesMostGapsToQuery(currentAlignment, '6QA2_A', env_variables, homologousSequences)
   breakpoint()
   filtered = list(filter(lambda al: al.seq == sequence_provides_most_gaps.seq,aligment_filtered))
   notFiltered = list(filter(lambda al: al.seq == sequence_provides_most_gaps.seq,currentAlignment))
   
   sequences_without_seq_most_gaps = functions.filterSequenceWithMostGaps(aligment_filtered, query_sec_header, env_variables, homologousSequences)
   seqs_homologous_after_load= homologousInCollection(sequences_without_seq_most_gaps,homologousSequences)
   #ver que aca tendria que agregar la homologa .
   assert len(seqs_homologous_after_load)==0
   assert len(filtered) == 0 and len(notFiltered) == 1


def test_filter_sequence_most_gapped_inject_homologous():
   # se agrega una homologa porque llega al minimo de secuencias
   breakpoint()
   currentAlignment = functions.loadFile('alignment.fasta')
   query_sec_header = '6QA2_A'

   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,query_sec_header)
   env_variables = variables_service.getDictVariablesValues()
   variables_service.setVariableEnv(MIN_SEQUENCES,95)
   mi_seq = variables_service.getVariableIntEnv(MIN_SEQUENCES)
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables)
   seqs_homologous_before_load= homologousInCollection(currentAlignment,homologousSequences)
   breakpoint()
   assert len(seqs_homologous_before_load)==0
   #busco la que más gaps tiene
   sequence_provides_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_sec_header, env_variables)

   aligment_filtered = functions.filterSequenceThatProvidesMostGapsToQuery(currentAlignment, '6QA2_A', env_variables, homologousSequences)
   breakpoint()
   
   #sequences_without_seq_most_gaps = functions.filterSequenceWithMostGaps(aligment_filtered, query_sec_header, env_variables, homologousSequences)
   seqs_homologous_after_load= homologousInCollection(aligment_filtered,homologousSequences)
   #ver que aca tendria que agregar la homologa .
   assert len(seqs_homologous_after_load)==1



def homologousInCollection(sequenceCollection,homologousSeqCollection):
   seq_hom_found =set()
   for sqh in homologousSeqCollection:
         seq_hom_found = set((list(filter(lambda al: al.id == sqh.id,sequenceCollection))))+seq_hom_found
   return seq_hom_found