from rest_framework import serializers
from .models import Organism, RegulatoryInteraction, ProjectAnalysisRegistry, Profile, User

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

class CreateUserSerializer(serializers.Serializer):
    firstName = serializers.CharField(max_length=30,required=True, allow_blank=False, allow_null=False)
    lastName = serializers.CharField(max_length=150,required=True, allow_blank=False, allow_null=False)
    email = serializers.EmailField(required=True, allow_blank=False, allow_null=False)
    organization = serializers.CharField(max_length=200,required=True, allow_blank=False, allow_null=False)
    password = serializers.CharField(max_length=20,required=True, allow_blank=False, allow_null=False)
    country = serializers.CharField(max_length=100,required=True, allow_blank=False, allow_null=False)
    brazilianState = serializers.CharField(max_length=2,required=False, allow_blank=True, allow_null=True)

class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ['id', 'first_name', 'email', 'password']
        extra_kwargs = {'password':{'write_only': True}}

class TestingSerializer(serializers.Serializer):
    tf_locus_tag = serializers.ListField()
    tg_locus_tag = serializers.ListField()
    connections = serializers.ListField()
