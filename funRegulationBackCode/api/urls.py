from rest_framework.routers import DefaultRouter
from api.views import OrganismViewSet, RegulatoryInteractionViewSet, ProjectAnalysisRegistryViewSet
from api.views import UserViewSet, LoginViewSet, UserLoginViewSet, RefreshTokenViewSet, LogoutViewSet
from api.views import ConfirmAccountViewSet, PasswordTokenCheckViewSet, RequestPasswordResetEmail
from api.views import SetNewPasswordViewSet, TaskStatusViewSet, CalculateCentralityViewSet

app_name = 'api'

router = DefaultRouter(trailing_slash=False)
router.register(r'Organisms', OrganismViewSet)
router.register(r'RegulatoryInteraction', RegulatoryInteractionViewSet, basename='RegulatoryInteraction')
router.register(r'CreateGrn', ProjectAnalysisRegistryViewSet, basename='CreateGrn')
router.register(r'CalculateCentrality', CalculateCentralityViewSet, basename='CalculateCentrality')
router.register(r'^TaskStatus/$', TaskStatusViewSet, basename='TastStatus')
router.register(r'CreateProfile', UserViewSet, basename='Profile')
router.register(r'Login', LoginViewSet, basename='Login')
router.register(r'UserTest', UserLoginViewSet, basename='UserTest')
router.register(r'Refresh', RefreshTokenViewSet, basename='Refresh')
router.register(r'Logout', LogoutViewSet, basename='Logout')
router.register(r'ConfirmAccount', ConfirmAccountViewSet, basename='ConfirmAccount')
router.register(r'RequestResetPassword', RequestPasswordResetEmail, basename='RequestPasswordReset')
router.register(r'^PasswordReset/(?P<uidb64>.+)/(?P<token>.+)/', PasswordTokenCheckViewSet, basename='PasswordReset')
router.register(r'^PasswordResetComplete', SetNewPasswordViewSet, basename='PasswordResetComplete')

urlpatterns = router.urls