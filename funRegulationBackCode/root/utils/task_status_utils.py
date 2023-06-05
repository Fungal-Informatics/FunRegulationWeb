from django.db.models import Case, When, Q, BooleanField

from funRegulationTool import task_utils


def get_analysis_fields():
    """
    Returns the possible fields with automatic analysis implementation.
    """
    return ['rsat', 'proteinortho']


def get_feature_status_statistic_generic(obj, fields):
    """
    Helper for get_feature_status_statistic(obj).
    """
    for field in fields:
        obj = obj.annotate(**{
            'last_%s_analysis_in_progress' % field: Case(When(Q(
                Q(**{'last_%s_registry__isnull' % field: False}) &
                ~Q(**{'last_%s_registry__registry__task__status' % field: task_utils.get_task_status_error()}) &
                Q(
                    Q(**{'last_%s_registry__task_%s__isnull' % (field, field): True}) |
                    ~Q(**{'last_%s_registry__task_%s__status__in' % (field, field):
                              task_utils.get_task_status_finished()})
                )
            ), then=True), default=False, output_field=BooleanField())})
        obj = obj.annotate(**{'last_%s_analysis_ok' % field: Case(When(Q(
            Q(**{'last_%s_registry__isnull' % field: False}) &
            Q(**{'last_%s_registry__registry__task__status' % field: task_utils.get_task_status_success()}) &
            Q(**{'last_%s_registry__task_%s__isnull' % (field, field): False}) &
            Q(**{'last_%s_registry__task_%s__status' % (field, field): task_utils.get_task_status_success()}) &
            Q(**{'last_%s_registry__%s_analysed' % (field, field): True}) &
            Q(**{'last_%s_registry__%s_error__isnull' % (field, field): True})
        ), then=True), default=False, output_field=BooleanField())})
    return obj


def get_feature_status_statistic(obj):
    """
    Annotation useful for observing the last analysis for each property of a GeneFeature.
    """
    obj = get_feature_status_statistic_generic(obj, get_analysis_fields())
    return obj


def has_analysis_in_progress_filter(in_progress=True):
    """
    Filters a feature which has an analysis process in progress or not.
    """
    return Q(
                Q(last_rsat_analysis_in_progress=in_progress) |
                Q(last_proteinortho_analysis_in_progress=in_progress)
    )


def has_analysis_with_error_filter():
    """
    Filters a feature which has analysis errors.
    """
    return Q(
                Q(last_rsat_registry__isnull=False,
                  last_rsat_analysis_in_progress=False, last_rsat_analysis_ok=False) |
                Q(last_proteinortho_registry__isnull=False,
                  last_proteinortho_analysis_in_progress=False, last_proteinortho_analysis_ok=False)
    )
