from rest_framework import serializers
from .models import Organism, RegulatoryInteraction, ProjectAnalysisRegistry, Profile, User
from django.contrib.auth.tokens import PasswordResetTokenGenerator
from django.utils.http import urlsafe_base64_decode, urlsafe_base64_encode
from django.contrib.sites.shortcuts import get_current_site
from django.urls import reverse
from django.utils.encoding import smart_str,force_str, smart_bytes, DjangoUnicodeDecodeError
from rest_framework.exceptions import AuthenticationFailed


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
    
class ProjectAnalysisRegistrySerializer(serializers.Serializer):
    organism_accession = serializers.CharField(max_length=150, required=True, allow_blank=False, allow_null=False)
    rsat_analyse = serializers.BooleanField(required=True, allow_null=False)
    download_organism = serializers.BooleanField(required=True, allow_null=False)

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

class UniqueTFsSerializer(serializers.Serializer):
    tf_locus_tag = serializers.CharField(max_length=100)

class UniqueTGsSerializer(serializers.Serializer):
    tg_locus_tag = serializers.CharField(max_length=100)

class RequestPasswordEmailSerializer(serializers.Serializer):
    email = serializers.EmailField()

class SetNewPasswordSerializer(serializers.Serializer):
    password = serializers.CharField(min_length=6, write_only=True)
    token = serializers.CharField(min_length=1, write_only=True)
    uidb64 = serializers.CharField(min_length=1, write_only=True)

    class Meta:
        fields = ['password', 'token', 'uidb64']
    
    def validate(self, attrs):
        try:
            password = attrs.get('password')
            token = attrs.get('token')
            uidb64 = attrs.get('uidb64')

            id = force_str(urlsafe_base64_decode(uidb64))
            user = User.objects.get(id = id)

            if not PasswordResetTokenGenerator().check_token(user, token):
                raise AuthenticationFailed('The reset link is invalid', 401)
            
            user.set_password(password)
            user.save()
        except Exception as e:
            raise AuthenticationFailed('The reset link is invalid', 401)
        return super().validate(attrs)
