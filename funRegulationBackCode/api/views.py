from django.shortcuts import render
from api.serializers import OrganismSerializer
from rest_framework import viewsets, permissions, mixins
#from root.modelsLenz.models import Organism
from .models import Organism

class OrganismViewSet(mixins.ListModelMixin,viewsets.GenericViewSet):
    queryset = Organism.objects.all()[:8]
    serializer_class = OrganismSerializer
    permission_classes = [permissions.IsAuthenticated]
