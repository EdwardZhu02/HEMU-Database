from django.urls import path, include, re_path
from Mainapp import views

urlpatterns = [
    # ==============ASYNC FRAMEWORK==============
    # Testing task handler
    path('asynctasks/test', views.async_test, name='async_test'),
    # Socket for querying task status
    path('asynctasks/status', views.get_async_task_progress, name='async_progress_query'),
    # Async task manager
    path('asynctasks/onhold', views.render_asyc_onhold, name='async_onhold'),

    # ==============HOME PAGE AND STATIC PAGES==============
    # Home Page
    path('', views.init_scr),
    path('home', views.init_scr, name='home_page'),
    # Site manager
    path('home/site_manage', views.site_manage, name='site_manage'),
    # More
    path('auxdata', views.aux_data_render, name='aux_data'),
    # path('auxdata/links', views.aux_data_render, name='external_links'),

    # ==============MODULE1: GENOME AND GENE INFORMATION==============
    path('genome/structure', views.genome_gene_sturcture_shinyapprender, name='genome_gene_structure'),
    path('genome/function', views.genome_gene_function_query_async, name='genome_gene_function'),
    path('genome/synteny', views.genome_gene_synteny_render, name='genome_gene_synteny'),
    path('databrowse', views.jbrowse_catalog_render, name='jbrowse'),
    path('blast', views.seqserver, name='sequenceserver'),

    # ==============MODULE2: TRANSCRIPTOME-DERIVED ANALYSIS==============
    path('gene/exp', views.gene_expression_query_async, name='gene_expression'),
    path('gene/sequence/', views.gene_sequence_obtain_async, name='gene_sequence'),
    path('gene/dge', views.gene_differential_analysis_async, name='gene_DE'),
    # DGE static file display
    re_path(r'^gene/static/(?P<identifier_name>\w+)/(?P<file_name>\w+.\w+);ht=(?P<frame_height>\d+);wid=('
            r'?P<frame_width>\d+)$', views.load_DE_staticfile),
    path('gene/gokegg', views.gene_gokegg_enrichment_async, name='gene_gokegg'),
    path('gene/wgcna', views.gene_wgcna_analysis_async, name='gene_wgcna'),

    # ============MODULE3: GENE FAMILY ANALYSIS=============
    path('genefam/hmm', views.genefam_identification_hmm_async, name='genefam_hmm'),
    path('genefam/blastp', views.genefam_identification_blastp_async, name='genefam_blastp'),
    path('genefam/phylo', views.genefam_phylogenetic_analysis_async, name='genefam_phylo'),
    path('genefam/htmap', views.genefam_expheatmapgen_async, name='genefam_expheatmap'),

    # ==============MODULE4: TE-DERIVED ANALYSIS==============
    path('te/exp', views.te_expression_query_async, name='te_expression'),
    path('te/insdensity', views.te_insertion_shinyapprender, name='te_insertion'),
    path('te/inslocation', views.te_location_search_async, name='te_location'),

    # ==============MODULE5: EPIGENOME ANALYSIS==============
    path('epigenome/chip', views.chip_peakanno_analysis_async, name='chip_peakanno'),
    path('epigenome/atacpeak', views.atac_peak_location_search_async, name='atac_peakloc_search'),

]
