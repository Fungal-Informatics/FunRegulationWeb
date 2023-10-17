import logging
from smtplib import SMTPException
from celery import shared_task
from django.utils import timezone
from django_celery_results.models import TaskResult
from funRegulationTool.task_utils import FunRegulationBaseTask
from django.core.mail import EmailMessage

logger = logging.getLogger('main')

def send_email(data):
    send_email_task.apply_async(args=[data])

@shared_task(bind=True, name='send_email', base=FunRegulationBaseTask)
def send_email_task(self, data):
    try:
        email = EmailMessage(subject=data['email_subject'], body=data['email_body'], to=[data['to_email']])
        email.send()
    except SMTPException as e:
        logger.error('error sending email of registry')
        self.retry(exc=e, countdown=600, max_retries=2)