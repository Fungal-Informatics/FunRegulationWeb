"""
    Construct GRN Orthology
"""
def construct_grn_orthology(model_organism_id, target_organism_id):
    
    model_regulatory_interactions = select_model_regulatory_by_organism_id(model_organism_id)
    
    for model_regulatory in model_regulatory_interactions:
        tf_orthologs = select_orthologs_by_target_organism(model_regulatory.tf_locus_tag, target_organism_id)
        tg_orthologs = select_orthologs_by_target_organism(model_regulatory.tg_locus_tag, target_organism_id)
        if len(tf_orthologs)!=0 and len(tg_orthologs)!=0:
            for ortholog_tf in tf_orthologs:
                for ortholog_tg in tg_orthologs:
                    regulatory_interaction = RegulatoryInteraction(0,ortholog_tf.target_protein.locus_tag, ortholog_tg.target_protein.locus_tag, model_regulatory.regulatory_function,model_regulatory.pubmedid)
                    
                    insert_regulatory_interaction(regulatory_interaction)
                    
                    #update gene as TF
                    tf = select_gene_by_locus_tag(ortholog_tf.target_protein.locus_tag)
                    tf.is_tf = 'True'
                    update_gene(tf)