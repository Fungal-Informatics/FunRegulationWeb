from django.shortcuts import render
from api.serializers import OrganismSerializer, RegulatoryInteractionSerializer, ProjectAnalysisRegistrySerializer
from api.serializers import CreateUserSerializer, TestingSerializer, UserSerializer
from rest_framework import viewsets, permissions, mixins
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.authentication import get_authorization_header
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import AllowAny
from rest_framework.exceptions import APIException, AuthenticationFailed
from django.db import DatabaseError, transaction
from .models import Organism, RegulatoryInteraction, Profile, Gene
from django.contrib.auth.models import User
from .authentication import create_access_token, refresh_access_token, decode_access_token, decode_refresh_token

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
        serializer = CreateUserSerializer(data=request.data)
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
    
class LoginViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def post(self, request):
        user = User.objects.filter(email=request.data['email']).first()

        if not user:
            raise APIException('Invalid credentials!')
        
        if not user.check_password(request.data['password']):
            raise APIException('Invalid credentials!')
        
        access_token = create_access_token(user.id)
        refresh_token = refresh_access_token(user.id)

        response = Response()
        response.set_cookie(key='refreshToken', value=refresh_token, httponly=True)
        response.data = {'token': access_token}
        return response
    
class UserLoginViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def list(self, request):
        auth = get_authorization_header(request).split()

        if auth and len(auth) == 2:
            token = auth[1].decode('utf-8')
            id = decode_access_token(token)
            user = User.objects.filter(pk=id).first()

            return Response(UserSerializer(user).data)
        return AuthenticationFailed('unaunthenticated')

class RefreshTokenViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def post(self, request):
        refresh_token = request.COOKIES.get('refreshToken')
        id = decode_refresh_token(refresh_token)
        access_token = create_access_token(id)
        return Response({
            'token': access_token
        })
    
class LogoutViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def post(self, request):
        respose = Response()
        respose.delete_cookie(key="refreshToken")
        respose.data = {
            'message:' : 'success'
        }
        return respose