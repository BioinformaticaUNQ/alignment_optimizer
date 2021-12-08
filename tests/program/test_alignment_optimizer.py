import pytest
from project_alignment_optimizer.program import variables_service
from project_alignment_optimizer.program import functions


#def test_filtrado():

 #   alignment = ["AA--C","A---C","AABDC","---DC","-ABD-"]
    # querySequence = 0 por defecto, es decir "AA--C"
#    query_sequence= ["AA--C","A---C","AABDC","-ABD-"]
   # var_services = variables_service.getDictVariablesValues()

  #  response = functions.filterAlignment(alignment,query_sequence,var_services)
    # Saca la secuencia "---DC" por ser la de mayor cantidad de gaps y menor score respecto de la query
 #   assert response == ["AA--C","A---C","AABDC","-ABD-"]

#    response = functions.filterAlignment(alignment)
    # Saca la secuencia "A---C" por ser la de mayor cantidad de gaps
 #   assert response == ["AA--C","AABDC","-ABD-"]

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
   
   #busco cual es la secuencia con mayor candidad de gaps
   query_seq = functions.find_alignment_by_header(currentAlignment,'6QA2_A')
   env_variables = variables_service.getDictVariablesValues()
   homologousSequences = functions.getHomologousSequences(query_seq, currentAlignment, env_variables)
   #valido que esa despues del primer filtrado no se encuentra 
   aligment_filtered = functions.filterSequenceThatProvidesMostGapsToQuery(currentAlignment, '6QA2_A', env_variables, homologousSequences)
   
  