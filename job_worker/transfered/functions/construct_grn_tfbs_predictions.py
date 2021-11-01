"""
    Construct GRN TFBS Predictions
"""
def construct_grn_tfbs_predictions(organism_id):
    tf_list = list()
    pwm_list = list()
    regulatory_interactions_list = list()
    
    # parse fasta file and turn into dictionary
    genome = SeqIO.to_dict(SeqIO.parse(open(in_file_genome), 'fasta'))
    
    # select TFs
    tf_list = select_tfs_by_organism(organism_id)
    for tf in tf_list:
        # Find PWMs
        pwm_list = select_pwms_by_locus_tag(tf.locus_tag)
        for pwm in pwm_list:
            #Find Regulatory Interactions
            regulatory_interactions_list = select_regulatory_interactions_by_tf_locus_tag(pwm.locus_tag)
            for regulatory_interaction in regulatory_interactions_list:
                
                prediction_file = out_folder + "tfbs_predictions/" + organism_id + "/" +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt"
                # this condition verify if exists a TFBS prediction already performed for this pwm
                # to avoid duplicated predictions with the same tg_promoter and pwm
                # note that pwms in Cis-Bp could be inferred by similarity tf_status = I
                if not os.path.isfile(prediction_file):
                    #Extract Promoter
                    promoter_sequence = extract_promoter(genome, regulatory_interaction.tg_locus_tag)
                    
                    ## Prepare matrix
                    in_matrix = in_folder + "Cis-BP/pwms/" + pwm.motif_id + ".txt"

                    # Input matrix (file content)
                    in_matrix_file = open(in_matrix, "r")
                    matrix = in_matrix_file.read()
                    in_matrix_file.close()

                    #Perform TFBS Prediction
                    ## Perform SOAP request on RSAT server with the matrix and the FASTA sequence
                    lib.log.info("Call RSAT Web Service: "+'\n'+
                             " tf_locus_tag: "+ regulatory_interaction.tf_locus_tag+'\n'+
                             " tg_locus_tag: "+ regulatory_interaction.tg_locus_tag+'\n'+
                             " pwm_id: "+ pwm.motif_id)

                    # Call web service matrix_scan
                    if matrix == ("Pos	A	C	G	T"+'\n'):
                        lib.log.info("Null Matrix: " + pwm.motif_id)
                    else:
                        result = call_matrix_scan(rsat_service, promoter_sequence, matrix)
                        # Write result in output file
                        with open(prediction_file, 'w') as out_file:
                            out_file.write(result)
                        out_file.close()

                        with open (prediction_file) as in_file:
                            for line in in_file:
                                if line.startswith("#"): continue
                                parts = line.strip().split("\t")
                                #Create new TFBS prediction for each RSAT prediction result
                                strand = urllib.parse.unquote(parts[3])
                                start = urllib.parse.unquote(parts[4])
                                end = urllib.parse.unquote(parts[5])
                                sequence = urllib.parse.unquote(parts[6])
                                weight = urllib.parse.unquote(parts[7])
                                pval = urllib.parse.unquote(parts[8])
                                ln_pval = urllib.parse.unquote(parts[9])
                                sig = urllib.parse.unquote(parts[10])
                                tfbs = TFBS(0, regulatory_interaction.id, pwm.id, strand, start, end, sequence, weight, pval, ln_pval, sig)
                                insert_tfbs_prediction(tfbs)
                            in_file.close()
                else:
                    lib.log.info("TFBS Prediction File already exists: " +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt")
