from rest_framework.routers import DefaultRouter
from api.views import OrganismViewSet, RegulatoryInteractionViewSet, ProjectAnalysisRegistryViewSet
from api.views import UserViewSet

app_name = 'api'

router = DefaultRouter(trailing_slash=False)
router.register(r'Organisms', OrganismViewSet)
router.register(r'^RegulatoryInteraction/$', RegulatoryInteractionViewSet, basename='RegulatoryInteraction')
router.register(r'ProjectAnalysisRegistry', ProjectAnalysisRegistryViewSet, basename='ProjectAnalysisRegistry')
router.register(r'CreateProfile', UserViewSet, basename='Profile')

urlpatterns = router.urls