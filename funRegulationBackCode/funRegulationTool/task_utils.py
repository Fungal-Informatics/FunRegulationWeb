import celery
from django_celery_results.models import TaskResult
from celery.states import STARTED

class FunRegulationBaseTask(celery.Task):
    def apply_async(self, *args, **kwargs):
        result = super(FunRegulationBaseTask, self).apply_async(*args, **kwargs)
        task_result, created = TaskResult.objects.get_or_create(task_id=result.task_id)
        task_result.status = result.state
        task_result.save()
        return result

    def __call__(self, *args, **kwargs):
        task_result, created = TaskResult.objects.get_or_create(task_id = self.request.id)
        task_result.status = STARTED
        task_result.save()
        return celery.Task.__call__(self, *args, **kwargs)


def get_task_status_finished():
    return ['SUCCESS', 'FAILURE']


def get_task_status_success():
    return 'SUCCESS'


def get_task_status_error():
    return 'FAILURE'