import logging

from celery import shared_task
from django_celery_results.models import TaskResult

from funRegulationTool.task_utils import FunRegulationBaseTask
from root.engine.calculate_centrality import *
from api.models import ProjectAnalysisRegistry, RegulatoryInteraction, CalculateCentralityRegistry, Gene, CalculateCentralityErrorType
import json
from networkx.readwrite import json_graph

logger = logging.getLogger('main')

@shared_task(bind=True, name='run_create_graph', base=FunRegulationBaseTask)
def task_run_create_graph(self, registry_id, organism_accession):
    logger.info(f'Running create graph for organism {organism_accession}')
    tfs = list()
    tgs = list()
    try:
        registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
        registry.task_create_graph = TaskResult.objects.get(task_id=self.request.id)
        registry.save()
    
        all_nodes = RegulatoryInteraction.objects.filter(organism_accession=organism_accession)
        
        for item in all_nodes:
            tfs.append(item.tf_locus_tag.locus_tag)
            tgs.append(item.tg_locus_tag.locus_tag)
        
        graph = create_graph(tf_nodes=tfs, tg_nodes=tgs)
        #ITS NECESSARY DO THIS TO PASS AS ARGUMENT IN CELERY 
        serialized_graph = json_graph.node_link_data(graph)
        logger.info(f'Create graph finished for organism {organism_accession}')
    except (Exception) as error:
        logger.error(f'An error ocurred in create graph task for {organism_accession}, error: {error}')
        registry.create_graph_error = CalculateCentralityErrorType.COMMAND_ERROR.value
        registry.save(update_fields=['create_graph_error'])
    return serialized_graph

@shared_task(bind=True, name='run_degree_centrality', base=FunRegulationBaseTask)
def task_degree_centrality(self, graph, registry_id, organism_accession):
    logger.info(f'Running calculate degree centrality for organism {organism_accession}')
    try:
        registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
        registry.task_calculate_degree = TaskResult.objects.get(task_id=self.request.id)
        registry.save()
        deserialized_graph = json_graph.node_link_graph(graph)
        results = calculate_degree_centrality(deserialized_graph)
        for gene_locus_tag, value in results.items():        
            centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                    degree = value)
            centrality.save()
        logger.info(f'End of calculate degree centrality for organism {organism_accession}')
    except (Exception) as error:
        logger.error(f'An error ocurred in degree calculate task for {organism_accession}, error: {error}')
        registry.degree_error = CalculateCentralityErrorType.COMMAND_ERROR.value
        registry.save(update_fields=['degree_error'])
    return graph

@shared_task(bind=True, name='run_closeness_centrality', base=FunRegulationBaseTask)
def task_closeness_centrality(self, graph, registry_id, organism_accession):
    logger.info(f'Running calculate closeness centrality for organism {organism_accession}')
    try:
        registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
        registry.task_calculate_closeness = TaskResult.objects.get(task_id=self.request.id)
        registry.save()
        deserialized_graph = json_graph.node_link_graph(graph)
        results = calculate_closeness_centrality(deserialized_graph)
        for gene_locus_tag, value in results.items():
            centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                    closeness = value)
            centrality.save(update_fields=['closeness'])
        logger.info(f'End of calculate closeness centrality for organism {organism_accession}')
    except (Exception) as error:
        logger.error(f'An error ocurred in closeness calculate task for {organism_accession}, error: {error}')
        registry.closeness_error = CalculateCentralityErrorType.COMMAND_ERROR.value
        registry.save(update_fields=['closeness_error'])
    return graph

@shared_task(bind=True, name='run_betweenness_centrality', base=FunRegulationBaseTask)
def task_betweenness_centrality(self, graph, registry_id, organism_accession):
    logger.info(f'Running calculate betweenness centrality for organism {organism_accession}')
    try:
        registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
        registry.task_calculate_betweenness = TaskResult.objects.get(task_id=self.request.id)
        registry.save()
        deserialized_graph = json_graph.node_link_graph(graph)
        results = calculate_betweenness_centrality(deserialized_graph)
        for gene_locus_tag, value in results.items():
            centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                    betweenness = value)
            centrality.save(update_fields=['betweenness'])
        logger.info(f'End of calculate betweenness centrality for organism {organism_accession}')
    except (Exception) as error:
        logger.error(f'An error ocurred in betweenness calculate task for {organism_accession}, error: {error}')
        registry.betweenness_error = CalculateCentralityErrorType.COMMAND_ERROR.value
        registry.save(update_fields=['betweenness_error'])
    return graph

@shared_task(bind=True, name='run_eigenvector_centrality', base=FunRegulationBaseTask)
def task_eigenvector_centrality(self, graph, registry_id, organism_accession):
    logger.info(f'Running calculate eigenvector centrality for organism {organism_accession}')
    try:
        registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
        registry.task_calculate_eigenvector = TaskResult.objects.get(task_id=self.request.id)
        registry.save()
        deserialized_graph = json_graph.node_link_graph(graph)
        results = calculate_eigenvector_centrality(deserialized_graph)
        for gene_locus_tag, value in results.items():
            centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                    eigenvector = value)
            centrality.save(update_fields=['eigenvector'])
        logger.info(f'End of calculate eigenvector centrality for organism {organism_accession}')
    except (Exception) as error:
        logger.error(f'An error ocurred in eigenvector calculate task for {organism_accession}, error: {error}')
        registry.eigenvector_error = CalculateCentralityErrorType.COMMAND_ERROR.value
        registry.save(update_fields=['eigenvector_error'])
    return graph

@shared_task(bind=True, name='run_harmonic_centrality', base=FunRegulationBaseTask)
def task_harmonic_centrality(self, graph, registry_id, organism_accession):
    logger.info(f'Running calculate harmonic centrality for organism {organism_accession}')
    try:
        registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
        registry.task_calculate_harmonic = TaskResult.objects.get(task_id=self.request.id)
        registry.save()
        deserialized_graph = json_graph.node_link_graph(graph)
        results = calculate_harmonic_centrality(deserialized_graph)
        for gene_locus_tag, value in results.items():
            centrality = CalculateCentralityRegistry(locus_tag = Gene.objects.get(locus_tag = gene_locus_tag), 
                                                    harmonic = value)
            centrality.save(update_fields=['harmonic'])
        logger.info(f'End of calculate harmonic centrality for organism {organism_accession}')
    except (Exception) as error:
        logger.error(f'An error ocurred in harmonic calculate task for {organism_accession}, error: {error}')
        registry.harmonic_error = CalculateCentralityErrorType.COMMAND_ERROR.value
        registry.save(update_fields=['harmonic_error'])
    return graph