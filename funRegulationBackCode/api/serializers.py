from rest_framework import serializers
from .models import Organism, RegulatoryInteraction

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

class RegulatoryInteractionSerializer(serializers.ModelSerializer):

    class Meta:
        model = RegulatoryInteraction
        fields = [
            'tf_locus_tag',
            'tg_locus_tag',
            'regulatory_function',
            'pubmedid_source'
        ]