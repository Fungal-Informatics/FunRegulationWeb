from api.serializers import OrganismSerializer, RegulatoryInteractionSerializer, ProjectAnalysisRegistrySerializer
from api.serializers import CreateUserSerializer, UniqueTFsSerializer, UniqueTGsSerializer, UserSerializer
from api.serializers import RequestPasswordEmailSerializer, SetNewPasswordSerializer, SearchGrnSerializer
from api.serializers import ProteinSerializer, CentralitySerializer
from rest_framework import viewsets, permissions, mixins, generics
from rest_framework.response import Response
from rest_framework import status
from rest_framework.authentication import get_authorization_header
from rest_framework.exceptions import APIException, AuthenticationFailed
from django.db import DatabaseError, transaction
from .models import Organism, RegulatoryInteraction, Profile, Gene, ProjectAnalysisRegistry, Protein, CalculateCentralityRegistry
from django.contrib.auth.models import User
from rest_framework_simplejwt.tokens import RefreshToken
from root.utils.tasks_email import send_email
from root.engine.calculate_centrality import *
from projects import tasks_chain
from django.contrib.sites.shortcuts import get_current_site
from django.urls import reverse
from django.conf import settings
from django_celery_results.models import TaskResult
import jwt
from datetime import datetime
from smtplib import SMTPException
from .authentication import create_access_token, refresh_access_token, decode_access_token,decode_refresh_token
from django.utils.http import urlsafe_base64_decode, urlsafe_base64_encode
from django.contrib.auth.tokens import PasswordResetTokenGenerator
from django.utils.encoding import smart_str,force_str, smart_bytes, DjangoUnicodeDecodeError


class OrganismViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    queryset = Organism.objects.all().filter(is_model=False).order_by('accession')
    serializer_class = OrganismSerializer
    permission_classes = [permissions.IsAuthenticated]

class RegulatoryInteractionViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    permission_classes = [permissions.IsAuthenticated]
    def post(self, request, format=None):   
        serializer = SearchGrnSerializer(data=request.data)
        if(serializer.is_valid()):
            data = serializer.data
            accession = data['organism_accession']
            queryset = RegulatoryInteraction.objects.select_related('tf_locus_tag')\
                .filter(tf_locus_tag__organism_accession=accession)\
                .distinct()\
                .order_by('tf_locus_tag','tg_locus_tag')
            
            unique_tfs_query = RegulatoryInteraction.objects.select_related('tf_locus_tag')\
                .filter(tf_locus_tag__organism_accession=accession)\
                .distinct()\
                .values('tf_locus_tag')\
                .order_by('tf_locus_tag','tg_locus_tag').distinct('tf_locus_tag')
            
            unique_tgs_query = RegulatoryInteraction.objects.select_related('tf_locus_tag')\
                .filter(tf_locus_tag__organism_accession=accession)\
                .distinct()\
                .values('tg_locus_tag')\
                .order_by('tg_locus_tag','tf_locus_tag').distinct('tg_locus_tag')
            
            unique_tfs = UniqueTFsSerializer(unique_tfs_query, many=True)
            unique_tgs = UniqueTGsSerializer(unique_tgs_query, many=True)
            serializer = RegulatoryInteractionSerializer(queryset, many=True)
            
            return Response({'tfs': unique_tfs.data, 'tgs': unique_tgs.data, 'edges':serializer.data}, status=status.HTTP_200_OK)

class ProjectAnalysisRegistryViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def post(self, request, format=None):
        serializer = ProjectAnalysisRegistrySerializer(data=request.data)
        if(serializer.is_valid()):
            with transaction.atomic():
                data = serializer.data
                registry = ProjectAnalysisRegistry(created_by=self.request.user)
                registry.organism_accession = data['organism_accession']
                registry.proteinortho_analyse = data['proteinOrtho_analyse']
                registry.rsat_analyse = data['rsat_analyse']
                registry.save()
                if((data['organism_accession'] == '') or (data['organism_accession'] is None)):
                    # UPLOAD FILES
                    registry.download_organism = False
                else:
                    registry.download_organism = True
                    tasks_chain.analyse_registry(registry)

            return Response({'registry': registry.pk, 'request': serializer.data}, status=status.HTTP_200_OK)
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
                    user.save()
                    profile = Profile.objects.create(user=user)
                    profile.organization = data['organization']
                    profile.BrazilianState = data['brazilianState']
                    profile.country = data['country']
                    profile.save()

                    # SEND EMAIL TO CONFIRM ACCOUNT
                    created_user = User.objects.get(email=data['email'])

                    token = RefreshToken.for_user(created_user).access_token

                    current_site = get_current_site(request).domain
                    #relative_link = reverse('ConfirmAccount-list')
                    relative_link = '/api/v1/ConfirmAccount'
                    absUrl = 'http://'+current_site+relative_link+"?token="+str(token)
                    email_body = 'Welcome '+ data['firstName']+ ' Use the link below to verify your account \n' + absUrl
                    email_data={'email_body':email_body, 'to_email':data['email'],
                                'email_subject': 'Verifiy your email'}
                    try:
                        send_email(email_data)
                    except SMTPException as e:
                        print(e)
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
    
class ConfirmAccountViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def list(self, request):
        ConfirmationToken = request.GET.get('token')
        try:
            payload = jwt.decode(ConfirmationToken, settings.SECRET_KEY, algorithms='HS256')
            user = Profile.objects.get(user=payload['user_id'])
            user.account_confirmation = True
            user.account_confirmation_date = datetime.now()
            user.save()
            return Response({'email': 'Success'}, status=status.HTTP_200_OK)
        except jwt.ExpiredSignatureError as e:
            return Response({'email': 'Activation Expired'}, status=status.HTTP_400_BAD_REQUEST)
        except jwt.exceptions.DecodeError as e:
            return Response({'email': 'Invalid Token'}, status=status.HTTP_400_BAD_REQUEST)

class RequestPasswordResetEmail(mixins.ListModelMixin, viewsets.GenericViewSet):
    serializer_class = RequestPasswordEmailSerializer

    def post(self, request):
        data = {'request': request,'data': request.data}
        serializer = self.serializer_class(data=data)
        email = request.data['email']
        if User.objects.filter(email = email).exists():
            user = User.objects.get(email = email)
            uidb64 = urlsafe_base64_encode(smart_bytes(user.id))
            token = PasswordResetTokenGenerator().make_token(user)
            current_site = get_current_site(request=request).domain
            #relative_link = reverse('api:PasswordReset-list', kwargs={'uidb64': uidb64, 'token': token})
            relative_link = '/api/v1/PasswordReset/'+str(uidb64)+'/'+str(token)+'/'
            absUrl = 'http://'+current_site+relative_link
            email_body = 'Hello, use the link below to reset your password \n' + absUrl
            email_data={'email_body':email_body, 'to_email':'d7073151@gmail.com','email_subject': 'Reset your password'}
            try:
                send_email(email_data)
            except SMTPException as e:
                print(e)
        return Response({'Success':'We sent the link to you reset your password'}, status=status.HTTP_200_OK)

class PasswordTokenCheckViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    def list(self, request, uidb64, token):
        try:
            #print(uidb64)
            id = smart_str(urlsafe_base64_decode(uidb64))
            user = User.objects.get(id=id)

            if not PasswordResetTokenGenerator().check_token(user, token):
                return Response({'error':'Token is not valid anymore, please request a new one'}, status=status.HTTP_401_UNAUTHORIZED)

            return Response({'success':True, 'message':'Credentials valid', 'uidb64': uidb64, 'token':token}, status=status.HTTP_200_OK)

        except DjangoUnicodeDecodeError as e:
            return Response({'error':'Token is not valid anymore, please request a new one'}, status=status.HTTP_401_UNAUTHORIZED)

class SetNewPasswordViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    serializer_class = SetNewPasswordSerializer

    def post(self, request):
        serializer = self.serializer_class(data = request.data)

        serializer.is_valid(raise_exception=True)
        return Response({'Success': True, 'Message': 'Password reset success'},status=status.HTTP_200_OK)
    
class TaskStatusViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    permission_classes = [permissions.IsAuthenticated]
    
    def list(self, request):
        id = self.request.query_params.get('id')
        index_header = 0
        tasks_completed = 0
        tasks = list()

        item = ProjectAnalysisRegistry.objects.filter(pk=id)\
        .values('rsat_analyse','task_id','task_download_organism_id','task_proteinortho_id','task_rsat_id',
                'task_create_graph_id','task_calculate_degree_id','task_calculate_closeness_id',
                'task_calculate_betweenness_id','task_calculate_eigenvector_id','task_calculate_harmonic_id')\
        .first()

        task_analyse = TaskResult.objects.filter(pk=item['task_id']).values('status').first()
        task_download = TaskResult.objects.filter(pk=item['task_download_organism_id']).values('status').first()
        task_proteinOrtho = TaskResult.objects.filter(pk=item['task_proteinortho_id']).values('status').first()
        task_create_graph = TaskResult.objects.filter(pk=item['task_create_graph_id']).values('status').first()
        task_calc_degree = TaskResult.objects.filter(pk=item['task_calculate_degree_id']).values('status').first()
        task_calc_closeness = TaskResult.objects.filter(pk=item['task_calculate_closeness_id']).values('status').first()
        task_calc_betweenness = TaskResult.objects.filter(pk=item['task_calculate_betweenness_id']).values('status').first()
        task_calc_eigenvector = TaskResult.objects.filter(pk=item['task_calculate_eigenvector_id']).values('status').first()
        task_calc_harmonic = TaskResult.objects.filter(pk=item['task_calculate_harmonic_id']).values('status').first()

        tasks.append(task_analyse)
        tasks.append(task_download)
        tasks.append(task_proteinOrtho)
        tasks.append(task_create_graph)
        tasks.append(task_calc_degree)
        tasks.append(task_calc_closeness)
        tasks.append(task_calc_betweenness)
        tasks.append(task_calc_eigenvector)
        tasks.append(task_calc_harmonic)

        if(item['rsat_analyse']):
            headers = ['Order analyse','Download organism','ProteinOrtho','Creating graph', 
                       'Degree calculation', 'Closeness calculation', 'Betweenness calculation',
                       'Eigenvector calculation', 'Harmonic calculation','Rsat', 'Pipeline Completed']
        
            result = {'Order analyse': '','Download organism': '','ProteinOrtho': '','Creating graph': '',
                    'Degree calculation': '','Closeness calculation': '','Betweenness calculation': '',
                    'Eigenvector calculation': '','Harmonic calculation': '','Rsat': '', 'Pipeline Completed': False}
            
            task_rsat = TaskResult.objects.filter(pk=item['task_rsat_id']).values('status').first()
            tasks.append(task_rsat)
        else:
            headers = ['Order analyse','Download organism','ProteinOrtho','Creating graph', 
                       'Degree calculation', 'Closeness calculation', 'Betweenness calculation',
                       'Eigenvector calculation', 'Harmonic calculation', 'Pipeline Completed']
            
            result = {'Order analyse': '','Download organism': '','ProteinOrtho': '','Creating graph': '',
                    'Degree calculation': '','Closeness calculation': '','Betweenness calculation': '',
                    'Eigenvector calculation': '','Harmonic calculation': '','Pipeline Completed': False}
        
        all_tasks = len(headers) - 1
        
        for task in tasks:
            if task is not None:
                if task['status'] == 'PENDING':
                    result[headers[index_header]] = 'WAITING FOR EXECUTION'
                elif task['status'] == 'STARTED':
                    result[headers[index_header]] = 'EXECUTING'
                elif task['status'] == 'FAILURE':
                    result[headers[index_header]] = 'ERROR'
                elif task['status'] == 'SUCCESS':
                    result[headers[index_header]] = 'COMPLETED'
                    tasks_completed += 1
            index_header += 1
        
        if(all_tasks == tasks_completed):
            result['Pipeline Completed'] = True if all_tasks == tasks_completed else False

        return Response(result, status=status.HTTP_200_OK)

class CalculateCentralityViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    permission_classes = [permissions.IsAuthenticated]
    def list(self, request):
        accession = self.request.query_params.get('accession')
        tfs = list()
        tgs = list()
        queryset = RegulatoryInteraction.objects.filter(organism_accession=accession)
        
        for item in queryset:
            tfs.append(item.tf_locus_tag.locus_tag)
            tgs.append(item.tg_locus_tag.locus_tag)
        
        graph = create_graph(tf_nodes=tfs, tg_nodes=tgs)
        degree = calculate_degree_centrality(graph)
        
        closeness = calculate_closeness_centrality(graph)
        betweenness = calculate_betweenness_centrality(graph)
        eigenvector = calculate_eigenvector_centrality(graph)
        harmonic = calculate_harmonic_centrality(graph)
        return Response({'degree': degree, 'closeness': closeness, 'betweenness': betweenness, 
                         'eigenvector': eigenvector, 'harmonic': harmonic}, status=status.HTTP_200_OK)
    
class ProteinInformationViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    permission_classes = [permissions.IsAuthenticated]
    def list(self, request):
        locus_tag = self.request.query_params.get('locus_tag')
        proteinInformation = Protein.objects.filter(locus_tag__locus_tag=locus_tag)
        serializer = ProteinSerializer(proteinInformation, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

class CentralityInformationViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    permission_classes = [permissions.IsAuthenticated]
    def list(self, request):
        locus_tag = self.request.query_params.get('locus_tag')
        centralityInformation = CalculateCentralityRegistry.objects.filter(locus_tag__locus_tag=locus_tag)
        serializer = CentralitySerializer(centralityInformation, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)
