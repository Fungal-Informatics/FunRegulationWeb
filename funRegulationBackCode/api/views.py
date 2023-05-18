from django.shortcuts import render
from api.serializers import OrganismSerializer, RegulatoryInteractionSerializer
from rest_framework import viewsets, permissions, mixins
from .models import Organism, RegulatoryInteraction

class OrganismViewSet(mixins.ListModelMixin,viewsets.GenericViewSet):
    queryset = Organism.objects.all()[:8]
    serializer_class = OrganismSerializer
    permission_classes = [permissions.AllowAny]

class RegulatoryInteractionViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    queryset = RegulatoryInteraction.objects.all()[:8]
    serializer_class = RegulatoryInteractionSerializer
    permission_classes = [permissions.AllowAny]
