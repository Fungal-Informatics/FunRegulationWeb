import logging

from celery import shared_task, group
from django.conf import settings
from django.db import transaction
from django_celery_results.models import TaskResult

from funRegulationTool.task_utils import FunRegulationBaseTask
from root.engine.calculate_centrality import *
from api.models import ProjectAnalysisRegistry, RegulatoryInteraction, CalculateCentralityRegistry, Gene
import json
from networkx.readwrite import json_graph
from time import sleep

@shared_task(bind=True, name='run_create_graph', base=FunRegulationBaseTask)
def task_run_create_graph(self, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_create_graph = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running create graph for organism %s' % organism_accession)
    tfs = list()
    tgs = list()

    all_nodes = RegulatoryInteraction.objects.filter(organism_accession=organism_accession)
        
    for item in all_nodes:
        tfs.append(item.tf_locus_tag.locus_tag)
        tgs.append(item.tg_locus_tag.locus_tag)
        
    graph = create_graph(tf_nodes=tfs, tg_nodes=tgs)
    logging.info('Create graph finished for organism %s' % organism_accession)
    serialized_graph = json_graph.node_link_data(graph)
    return serialized_graph

@shared_task(bind=True, name='run_degree_centrality', base=FunRegulationBaseTask)
def task_degree_centrality(self, graph, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_calculate_degree = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running calculate degree centrality for organism')
    graph1 = json_graph.node_link_graph(graph)
    results = calculate_degree_centrality(graph1)
    logging.info('End of calculate degree centrality for organism')
    for gene_locus_tag, value in results.items():        
        centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                 degree = value)
        centrality.save()
    logging.info('End of calculate degree centrality for organism %s' % organism_accession)
    return graph

@shared_task(bind=True, name='run_closeness_centrality', base=FunRegulationBaseTask)
def task_closeness_centrality(self, graph, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_calculate_closeness = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running calculate closeness centrality for organism %s' % organism_accession)
    graph1 = json_graph.node_link_graph(graph)
    results = calculate_closeness_centrality(graph1)
    logging.info('End of calculate closeness centrality for organism %s' % organism_accession)
    for gene_locus_tag, value in results.items():
        centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                 closeness = value)
        centrality.save(update_fields=['closeness'])
    logging.info('End of calculate closeness centrality for organism %s' % organism_accession)
    return graph

@shared_task(bind=True, name='run_betweenness_centrality', base=FunRegulationBaseTask)
def task_betweenness_centrality(self, graph, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_calculate_betweenness = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running calculate betweenness centrality for organism %s' % organism_accession)
    graph1 = json_graph.node_link_graph(graph)
    results = calculate_betweenness_centrality(graph1)
    for gene_locus_tag, value in results.items():
        centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                betweenness = value)
        centrality.save(update_fields=['betweenness'])
    logging.info('End of calculate betweenness centrality for organism %s' % organism_accession)
    return graph

@shared_task(bind=True, name='run_eigenvector_centrality', base=FunRegulationBaseTask)
def task_eigenvector_centrality(self, graph, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_calculate_eigenvector = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running calculate eigenvector centrality for organism %s' % organism_accession)
    graph1 = json_graph.node_link_graph(graph)
    results = calculate_eigenvector_centrality(graph1)
    for gene_locus_tag, value in results.items():
        centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                 eigenvector = value)
        centrality.save(update_fields=['eigenvector'])
    logging.info('End of calculate eigenvector centrality for organism %s' % organism_accession)
    return graph

@shared_task(bind=True, name='run_harmonic_centrality', base=FunRegulationBaseTask)
def task_harmonic_centrality(self, graph, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_calculate_harmonic = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running calculate harmonic centrality for organism %s' % organism_accession)
    graph1 = json_graph.node_link_graph(graph)
    results = calculate_harmonic_centrality(graph1)
    for gene_locus_tag, value in results.items():
        centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                 harmonic = value)
        centrality.save(update_fields=['harmonic'])
    logging.info('End of calculate harmonic centrality for organism %s' % organism_accession)
    return graph