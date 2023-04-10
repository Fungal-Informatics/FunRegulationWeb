from rest_framework.routers import DefaultRouter
from api.views import OrganismViewSet

app_name = 'api'

router = DefaultRouter(trailing_slash=False)
router.register(r'Organisms', OrganismViewSet)

urlpatterns = router.urls