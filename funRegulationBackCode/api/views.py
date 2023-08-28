from api.serializers import OrganismSerializer, RegulatoryInteractionSerializer, ProjectAnalysisRegistrySerializer
from api.serializers import CreateUserSerializer, UniqueTFsSerializer, UniqueTGsSerializer, UserSerializer
from api.serializers import RequestPasswordEmailSerializer, SetNewPasswordSerializer
from rest_framework import viewsets, permissions, mixins, generics
from rest_framework.response import Response
from rest_framework import status
from rest_framework.authentication import get_authorization_header
from rest_framework.exceptions import APIException, AuthenticationFailed
from django.db import DatabaseError, transaction
from .models import Organism, RegulatoryInteraction, Profile, Gene
from django.contrib.auth.models import User
from rest_framework_simplejwt.tokens import RefreshToken
from root.utils.tasks_email import send_email
from django.contrib.sites.shortcuts import get_current_site
from django.urls import reverse
from django.conf import settings
import jwt
from datetime import datetime
from smtplib import SMTPException
from .authentication import create_access_token, refresh_access_token, decode_access_token,decode_refresh_token
from django.utils.http import urlsafe_base64_decode, urlsafe_base64_encode
from django.contrib.auth.tokens import PasswordResetTokenGenerator
from django.utils.encoding import smart_str,force_str, smart_bytes, DjangoUnicodeDecodeError

class OrganismViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer
    permission_classes = [permissions.IsAuthenticated]

class RegulatoryInteractionViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    #permission_classes = [permissions.IsAuthenticated]
    permission_classes = [permissions.AllowAny]
    def list(self, request):        
        accession = self.request.query_params.get('organism_accession')
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
            print(uidb64)
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