import re

import unicodedata

from django.db.models import Case, When, BooleanField, Q


def get_str_without_formatting(str_original):
    """
    Gets string without special characters.
    Useful for creating file names.
    """
    if not str_original:
        return None
    str_original = str_original.strip()
    str_original = unicodedata.normalize('NFKD', str_original)
    str_original = u"".join([c for c in str_original if not unicodedata.combining(c)])
    str_original = str_original.lower()
    str_original = str_original.replace(' ', '_')
    str_original = ''.join(re.findall("[a-z0-9_.-]+", str_original))
    return str_original


def only_numbers(valor):
    """
    Gets only numbers from a string.
    """
    if valor and type(valor) is str:
        return ''.join(re.findall('\d+', valor))
    return None


def int_to_bool_annotation(field_name):
    return Case(When(Q(**{'%s__gte' % field_name: 1}), then=True), defaul=False, output_field=BooleanField())


def get_valid_feature_filters(project, prefix=''):
    return Q(
        ~Q(**{'%sfeature' % prefix: 'undefined'}),
        **{
            '%sremoved' % prefix: False,
            '%sgene__removed' % prefix: False,
            '%sgene__project' % prefix: project
        }
    )
