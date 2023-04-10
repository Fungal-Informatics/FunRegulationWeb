from rest_framework import serializers
#from root.modelsLenz.models import Organism
from .models import Organism

class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = [
            'accession',
            'order',
            'genus',
            'species',
            'strain',
            'is_model',
            'cis_bp'
        ]