from django.shortcuts import render
from api.serializers import OrganismSerializer, RegulatoryInteractionSerializer, ProjectAnalysisRegistrySerializer
from api.serializers import UserSerializer, TestingSerializer
from rest_framework import viewsets, permissions, mixins
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import AllowAny
from django.db import DatabaseError, transaction
from .models import Organism, RegulatoryInteraction, Profile, Gene
from django.contrib.auth.models import User
from itertools import chain

class OrganismViewSet(mixins.ListModelMixin,viewsets.GenericViewSet):
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer
    permission_classes = [permissions.IsAuthenticated]

class RegulatoryInteractionViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    permission_classes = [permissions.IsAuthenticated]
    def list(self, request):
        all_tf = list()
        all_tg = list()
        unique_tfs = list()
        unique_tgs = list()
        accession = self.request.query_params.get('organism_accession')
        queryset = RegulatoryInteraction.objects.select_related('tf_locus_tag')\
            .filter(tf_locus_tag__organism_accession=accession)\
            .distinct()\
            .order_by('tf_locus_tag','tg_locus_tag')
        for item in queryset:
            all_tf.append(item.tf_locus_tag.locus_tag)
            all_tg.append(item.tg_locus_tag.locus_tag)
        
        unique_tfs = [*set(all_tf)]
        unique_tgs = [*set(all_tg)]
        serializer = RegulatoryInteractionSerializer(queryset, many=True)
        return Response({'tfs': unique_tfs,'tgs': unique_tgs,'edges':serializer.data},status=status.HTTP_200_OK)


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
                if(not(User.objects.filter(username=self.request.user).count()>=1)):
                    data = serializer.data
                    user = User()
                    user.first_name = data['firstName']
                    user.last_name = data['lastName']
                    user.username = data['email']
                    user.set_password(data['password'])
                    user.email = user.username
                    #user.save()
                    profile = Profile.objects.create(user=user)
                    profile.organization = data['organization']
                    profile.BrazilianState = data['brazilianState']
                    profile.country = data['country']
                    #profile.save()
                    return Response(serializer.data,status=status.HTTP_201_CREATED)
                else:
                    data = serializer.data
                    user = User.objects.get(username=self.request.user)
                    user.first_name = data['firstName']
                    user.last_name = data['lastName']
                    user.username = data['email']
                    user.email = data['email']
                    user.set_password(data['password'])
                    user.save()
                    profile = Profile.objects.get(user_id=user.pk)
                    profile.organization = data['organization']
                    profile.BrazilianState = data['brazilianState']
                    profile.country = data['country']
                    profile.save()
                    return Response(serializer.data,status=status.HTTP_200_OK)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)