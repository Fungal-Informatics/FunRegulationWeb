"""
    Select orthologs by model_locus_tag and target organism id
"""
def select_orthologs_by_target_organism(model_locus_tag, target_organism_id):
    orthology_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT model_protein, target_protein from orthology ortho join protein prot on ortho.model_protein = prot.id join gene gen on prot.locus_tag = %s AND gen.organism = %s WHERE ortho.model_protein IS NOT NULL", (model_locus_tag, target_organism_id))
        rec = cursor.fetchall()
        for row in rec:
            ortho = Orthology(select_protein_by_id(row[0]),select_protein_by_id(row[1]))
            orthology_list.append(ortho)
        cursor.close()
        return orthology_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table orthology", error)
        lib.log.info(model_locus_tag)