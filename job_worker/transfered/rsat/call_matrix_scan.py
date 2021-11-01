"""
    RSAT Web Service: Scan sequences for a given matrix
"""
def call_matrix_scan(service, fasta_content_str, matrix_str):
    #RSAT config parameters
    uth_pval = '1e-2'
    matrix_format = "cis-bp"
    
    print(fasta_content_str)
    print(matrix_str)
    # Wrap all arguments into a named list
    #http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
    arguments = {
            'sequence' : fasta_content_str,
            'matrix' : matrix_str,
            'matrix_format' : matrix_format,
            'uth' : ['pval '+uth_pval],
            'lth' : ['score '+str(1)],
            #'quick' : 1,
            'str' : 2,
            'origin' : 'start',
            'background_input' : 1, # this option requires 'markov'
            'markov' : 1,
            'pseudo' : 1,
            'decimals' : 1,
            'n_treatment' : 'score',
            'background_pseudo': 0.01
            #'sequence_format' : 'fasta',
            #'organism': 'fungi'
            #'verbosity' : 1
    }

    ## Perform SOAP request on RSAT server
    result = service.matrix_scan(arguments)
    print (result.client)
    return result.client