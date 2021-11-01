"""
    Parse orthology file - ProteinOrtho
"""
def parse_orthology_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            line_parts = line.strip().split("\t")
            
            model = urllib.parse.unquote(line_parts[3])
            target = urllib.parse.unquote(line_parts[4])
            
            model_parts = model.strip().split(",")
            target_parts = target.strip().split(",")
            
            for record_model in model_parts:
                for record_target in target_parts:
                    if (record_model != '*' and record_target != '*'):
                        model_protein = select_protein_by_id(record_model)
                        target_protein = select_protein_by_id(record_target)
                        orthology = Orthology(model_protein,target_protein)
                        insert_orthology(orthology)
                    
    in_file.close()
    lib.log.info(filename + " parsed correctly")