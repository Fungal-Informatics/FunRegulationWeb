"""
    Parse Ensembl Species input file
"""
def parse_organism_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            order = urllib.parse.unquote(parts[0])
            genus = urllib.parse.unquote(parts[1])
            species = urllib.parse.unquote(parts[2])
            strain = urllib.parse.unquote(parts[3])
            taxon = urllib.parse.unquote(parts[4])
            assembly_name = urllib.parse.unquote(parts[5])
            assembly_accession = urllib.parse.unquote(parts[6])
            is_model = urllib.parse.unquote(parts[7])
            cis_bp = urllib.parse.unquote(parts[8])
            organism = Organism(order,genus,species,strain,taxon,assembly_name,assembly_accession,is_model,cis_bp)
            insert_organism(organism)
    in_file.close()
    lib.log.info(filename + " parsed correctly")