from rest_framework import serializers
from .models import Organism, RegulatoryInteraction, ProjectAnalysisRegistry, Profile

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
    
class ProjectAnalysisRegistrySerializer(serializers.ModelSerializer):
    class Meta:
        model = ProjectAnalysisRegistry
        fields = [
            'organism_accession',
            'rsat_analyse'
        ]

class UserSerializer(serializers.Serializer):
    firstName = serializers.CharField(max_length=30)
    lastName = serializers.CharField(max_length=150)
    email = serializers.EmailField()
    organization = serializers.CharField(max_length=200)
    password = serializers.CharField(max_length=20)
    country = serializers.CharField(max_length=100)
    brazilianState = serializers.CharField(max_length=2)
