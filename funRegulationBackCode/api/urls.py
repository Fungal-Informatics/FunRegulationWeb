from rest_framework.routers import DefaultRouter
from api.views import OrganismViewSet, RegulatoryInteractionViewSet

app_name = 'api'

router = DefaultRouter(trailing_slash=False)
router.register(r'Organisms', OrganismViewSet)
router.register(r'RegulatoryInteraction', RegulatoryInteractionViewSet)

urlpatterns = router.urls