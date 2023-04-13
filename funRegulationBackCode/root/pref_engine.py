import json
import logging

from django.core.exceptions import ObjectDoesNotExist
from django.db import DatabaseError
from django.db.models import Q

from api.models import SystemPreference, SystemPreferenceType


class PrefEngine:
    """
    Engine for system preferences.
    """
    @classmethod
    def get_pref(cls, name, default, user=None, project=None):
        """
        Gets a preference value.
        """
        try:
            filters = Q(name=name)
            if user:
                filters &= Q(user=user)
            if project:
                filters &= Q(project=project)
            #pref = SystemPreference.objects.get(filters)
            pref = None
            return cls.__convert_value_pref(pref)
        except ObjectDoesNotExist:
            logging.debug('preference %s for user %s and project %s does not exist. Applying default value %s' %
                          (name, user, project, default))
        except DatabaseError:
            logging.exception('error loading preference %s for user %s and project %s. Applying default value %s' %
                              (name, user, project, default))
        return default

    @classmethod
    def set_pref(cls, name, pref_type, value, user=None, project=None):
        """
        Sets a preference value.
        """
        try:
            pref = SystemPreference.objects.get(name=name, user=user, project=project)
        except ObjectDoesNotExist:
            pref = SystemPreference()
            pref.name = name
            pref.pref_type = pref_type
            pref.user = user
            pref.project = project
        pref.valor = str(value)
        pref.save()

    @classmethod
    def __convert_value_pref(cls, preference):
        """
        Converts a preference value to its real type.
        """
        if preference.pref_type == SystemPreferenceType.INTEGER.value:
            return int(preference.value)
        if preference.pref_type == SystemPreferenceType.FLOAT.value:
            return float(preference.value)
        if preference.pref_type == SystemPreferenceType.JSON.value:
            return json.loads(preference.value)
        if preference.pref_type == SystemPreferenceType.BOOLEAN.value:
            return preference.value.lower() == 'true'
        return preference.value
