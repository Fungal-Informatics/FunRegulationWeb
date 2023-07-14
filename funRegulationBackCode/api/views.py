from django.shortcuts import render
from api.serializers import OrganismSerializer, RegulatoryInteractionSerializer, ProjectAnalysisRegistrySerializer
from api.serializers import UserSerializer
from rest_framework import viewsets, permissions, mixins
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import AllowAny
from django.db import DatabaseError, transaction
from .models import Organism, RegulatoryInteraction, Profile
from django.contrib.auth.models import User

class OrganismViewSet(mixins.ListModelMixin,viewsets.GenericViewSet):
    queryset = Organism.objects.all()[:8]
    serializer_class = OrganismSerializer
    permission_classes = [permissions.AllowAny]

class RegulatoryInteractionViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    queryset = RegulatoryInteraction.objects.all()[:8]
    serializer_class = RegulatoryInteractionSerializer
    permission_classes = [permissions.AllowAny]

class ProjectAnalysisRegistryViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def post(self, request, format=None):
        serializer = ProjectAnalysisRegistrySerializer(data=request.data)
        if(serializer.is_valid()):
            # with transaction.atomic:
            # accession = self.validated_data.get("organism_accession")
            # rsat = self.validated_data.get("rsat_analyse")
            serializer.organism_accession = "ENTREI AQUI"
            #serializer.save()
            return Response(serializer.data,status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class UserViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def post(self, request, format=None):
        serializer = UserSerializer(data=request.data)
        if(serializer.is_valid()):
            with transaction.atomic():
                data = serializer.data
                user = User()
                user.first_name = data['firstName']
                user.last_name = data['lastName']
                user.username = data['email']
                user.set_password(data['password'])
                user.email = user.username
                user.save()
                profile = Profile.objects.create(user=user)
                profile.organization = data['organization']
                profile.BrazilianState = data['brazilianState']
                profile.country = data['country']
                profile.save()
                return Response(serializer.data,status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)