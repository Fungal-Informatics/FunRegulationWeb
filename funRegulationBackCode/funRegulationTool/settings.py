import logging
import os
from decouple import config
from django.core.files.storage import FileSystemStorage

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/2.1/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = ')ir)4m(0nhwwad*v-8^q3yf)yn2u#n%9@$hwi8a=w*3+j$4b!j'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ['*']


# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_celery_results',
    'root',
    'projects',
    'rest_framework',
    'api'
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'funRegulationTool.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, 'main/templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'funRegulationTool.wsgi.application'


# Database
# https://docs.djangoproject.com/en/2.1/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': config('FUNREGULATION_DATABASE_NAME'),
        'USER': config('FUNREGULATION_DATABASE_USER'),
        'PASSWORD': config('FUNREGULATION_DATABASE_PASSWORD'),
        'HOST': config('FUNREGULATION_DATABASE_HOST'),
        'PORT': config('FUNREGULATION_DATABASE_PORT')
    }
}


# Password validation
# https://docs.djangoproject.com/en/2.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/2.1/topics/i18n/

LANGUAGE_CODE = 'en-us'
TIME_ZONE = config('FUNREGULATION_TIME_ZONE')
USE_I18N = True
USE_L10N = True
USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.1/howto/static-files/

STATIC_URL = '/static/'

CELERY_BROKER_URL = 'pyamqp://guest@localhost//'
CELERY_RESULT_BACKEND = 'django-db'

CELERY_TASK_ROUTES = {
    #task : queue
    'import_genes': 'import_genes',
    'analyse_registry': 'analyse_registry',
    'run_proteinortho': 'run_proteinortho',
    'run_rsat': 'run_rsat',
}

# GENERAL FILES
PWMS_FILES_PATH = config('FUNREGULATION_PWMS_FILES')
LOG_FILE_PATH = config('FUNREGULATION_LOG_FILE')
NCBI_DOWNLOAD_PATH = config('FUNREGULATION_DOWNLOAD_NCBI_PATH')

# MODEL ORGANISMS FILES
ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PATH = config('FUNREGULATION_ORGANISM_MODEL_FUSARIUM_GRAMINEARUM')
ORGANISM_MODEL_NEUROSPORA_CRASSA_PATH = config('FUNREGULATION_ORGANISM_MODEL_NEUROSPORA_CRASSA')
ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PATH = config('FUNREGULATION_ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE')
ORGANISM_MODEL_A_NIDULANS_PATH = config('FUNREGULATION_ORGANISM_MODEL_A_NIDULANS')

# MODEL ORGANISMS PROTEIN FILES
ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PROTEIN_PATH = config('FUNREGULATION_ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PROTEIN')
ORGANISM_MODEL_NEUROSPORA_CRASSA_PROTEIN_PATH = config('FUNREGULATION_ORGANISM_MODEL_NEUROSPORA_CRASSA_PROTEIN')
ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PROTEIN_PATH = config('FUNREGULATION_ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PROTEIN')
ORGANISM_MODEL_A_NIDULANS_PROTEIN_PATH = config('FUNREGULATION_ORGANISM_MODEL_A_NIDULANS_PROTEIN')

# RSAT
RSAT_PATH = config('FUNREGULATION_RSAT_PATH')
RSAT_SAVE_PATH = config('FUNREGULATION_RSAT_SAVE_PATH')

# PROTEINORTHO
PROTEINORTHO_PATH = config('FUNREGULATION_PROTEINORTHO_PATH')
PROTEINORTHO_SAVE_PATH = config('FUNREGULATION_PROTEINORTHO_SAVE_PATH')