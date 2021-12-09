import pytest
from project_alignment_optimizer.program import variables_service
from project_alignment_optimizer.program import functions

def test_search_query_sequence():
    breakpoint()
    file = functions.loadFile('alignment.fasta')
    #print( file)
    
    #print('var_services ' + str(var_services))
    #variables_service.setVariableEnv(variables_service.QUERY_SEQUENCE_HEADER, '6QA2_A')
    #variables_service.dotenv.load_dotenv(variables_service.dotenv_file, override=True) # Need to reload de os.env variables
    breakpoint()
    print('el header ' + '6QA2_A')
    query_seq = functions.find_alignment_by_header(file,'6QA2_A')
    print('la query_seq ' + str(query_seq))
    ungapped = query_seq.seq.ungap()
    #print('las queries' + ungapped)
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
  
def test_filter_sequence_most_gapped_inject_homologus_found():
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
   sequence_with_most_gaps = functions.sequenceProvidesMostGaps(currentAlignment, query_sec_header, env_variables)
   #aca ver que traiga las seqs no la que va a sacar 
   breakpoint()
   aligment_filtered = functions.filterSequenceThatProvidesMostGapsToQuery(currentAlignment, '6QA2_A', env_variables, homologousSequences)
   breakpoint()
   filtered = list(filter(lambda al: al.seq == sequence_with_most_gaps.seq,aligment_filtered))
   notFiltered = list(filter(lambda al: al.seq == sequence_with_most_gaps.seq,currentAlignment))
   #ver cuando agrega homologas
   functions.filterSequenceWithMostGaps(aligment_filtered, query_sec_header, env_variables, homologousSequences)
   seqs_homologous_after_load= homologousInCollection(sequence_with_most_gaps,homologousSequences)
   #ver que aca tendria que agregar la homologa .
   assert len(seqs_homologous_after_load)==1
   assert len(filtered) == 0 and len(notFiltered) == 1




def homologousInCollection(sequenceCollection,homologousSeqCollection):
   seq_hom_found =[]
   for sqh in homologousSeqCollection:
         seq_hom_found+(list(filter(lambda al: al.seq.ungap() == sqh.seq,sequenceCollection)))
   return seq_hom_found