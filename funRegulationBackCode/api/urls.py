from rest_framework.routers import DefaultRouter
from api.views import OrganismViewSet, RegulatoryInteractionViewSet, ProjectAnalysisRegistryViewSet
from api.views import UserViewSet, LoginViewSet, UserLoginViewSet, RefreshTokenViewSet, LogoutViewSet

app_name = 'api'

router = DefaultRouter(trailing_slash=False)
router.register(r'Organisms', OrganismViewSet)
router.register(r'^RegulatoryInteraction/$', RegulatoryInteractionViewSet, basename='RegulatoryInteraction')
router.register(r'ProjectAnalysisRegistry', ProjectAnalysisRegistryViewSet, basename='ProjectAnalysisRegistry')
router.register(r'CreateProfile', UserViewSet, basename='Profile')
router.register(r'Login', LoginViewSet, basename='Login')
router.register(r'UserTest', UserLoginViewSet, basename='UserTest')
router.register(r'Refresh', RefreshTokenViewSet, basename='Refresh')
router.register(r'Logout', LogoutViewSet, basename='Logout')

urlpatterns = router.urls