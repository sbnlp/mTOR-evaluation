#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import datetime
import glob
import networkx
import numpy
import os
import pandas
import pickle

import networkx_analysis

nodes_overlap_species_nm_fns = [ networkx_analysis.nm_name_equal, 
                   networkx_analysis.nm_name_clean_equal, 
                   networkx_analysis.nm_name_clean_approx, 
                   #networkx_analysis.nm_gene_id_intersect, 
                   #networkx_analysis.nm_name_clean_approx_OR_gene_id_intersect,
                   networkx_analysis.nm_bqbiol_is_equal,
                   networkx_analysis.nm_bqbiol_is_overlaps,
                   
                   networkx_analysis.nm_name_equal_w_participants, 
                   networkx_analysis.nm_name_clean_equal_w_participants, 
                   networkx_analysis.nm_name_clean_approx_w_participants, 
                   #networkx_analysis.nm_gene_id_intersect_w_participants, 
                   networkx_analysis.nm_bqbiol_is_equal_w_participants,
                   networkx_analysis.nm_bqbiol_is_overlaps_w_participants]
                   
nodes_overlap_species_nm_fns_paper = [ 
                   networkx_analysis.nm_name_clean_equal, 
                   networkx_analysis.nm_name_clean_approx, 
                   networkx_analysis.nm_bqbiol_is_equal,
                   networkx_analysis.nm_bqbiol_is_overlaps,
                   networkx_analysis.nm_name_clean_equal_w_participants, 
                   networkx_analysis.nm_name_clean_approx_w_participants, 
                   networkx_analysis.nm_bqbiol_is_equal_w_participants,
                   networkx_analysis.nm_bqbiol_is_overlaps_w_participants]
nodes_overlap_species_nm_fns_paper_names = ["nmeq", "appeq", "enteq", "entov", "nmeq/wc", "appeq/wc", "enteq/wc", "entov/wc"]

nodes_overlap_reactions_nm_fns = [ networkx_analysis.nm_bqbiol_is_equal, 
                                  networkx_analysis.nm_bqbiol_is_overlaps, 
                                  networkx_analysis.nm_bqbiol_is_overlaps_sbo_is_a]
nodes_overlap_reactions_nm_fns_names = ["sboeq", "sboov", "sboisa"]
                                              
subgraphs_overlap_node_match_fns = [networkx_analysis.nm_name_clean_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                      networkx_analysis.nm_name_clean_approx_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_approx_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_approx_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_overlaps_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_overlaps_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_overlaps_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                   
                      networkx_analysis.nm_name_clean_w_participants_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_w_participants_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                      networkx_analysis.nm_name_clean_approx_w_participants_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_approx_w_participants_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_approx_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_equal_w_participants_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_equal_w_participants_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_equal_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_overlaps_w_participants_AND_nm_bqbiol_is_equal,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_overlaps_w_participants_AND_nm_bqbiol_is_overlaps,
                      networkx_analysis.nm_name_clean_approx_OR_bqbiol_is_overlaps_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a]  

subgraphs_overlap_node_match_fns_names = ["nmeq, sboeq", "nmeq, sboov", "nmeq, sboisa", 
                                          "appeq, sboeq", "appeq, sboov","appeq, sboisa",
                                          "appeq/enteq, sboeq", "appeq/enteq, sboov","appeq/enteq, sboisa",
                                          "appeq/entov, sboeq", "appeq/entov, sboov", "appeq/entov, sboisa",
                                          "nmeq/wc, sboeq", "nmeq/wc, sboov", "nmeq/wc, sboisa",
                                          "appeq/wc, sboeq", "appeq/wc, sboov", "appeq/wc, sboisa",
                                          "appeq/enteq/wc, sboeq", "appeq/enteq/wc, sboov", "appeq/enteq/wc, sboisa",
                                          "appeq/entov/wc, sboeq", "appeq/entov/wc, sboov", "appeq/entov/wc, sboisa"]
subgraphs_overlap_node_match_fns_names_map = { k.__name__: v for k, v in zip( subgraphs_overlap_node_match_fns, subgraphs_overlap_node_match_fns_names)}

########################################################################
########################################################################
# INITIALIZE

def mtor_dir():
    directory = os.environ.get("MTOR_DATA_DIR")
    if directory is None:
      directory = ""
    return directory + "/"

def mtor_dir_results():
    return mtor_dir() + "/results/"

def mtor_dir_results_graphs():
    return mtor_dir() + "/results-graphs/"

def mtor_dir_results_statistics():
    return mtor_dir() + "/results-statistics/"
    
def now():
    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

########################################################################
########################################################################
# HELPERS

def plot_precision_recall_f_score( data, name, x_ticks = None, columns = ["precision", "recall","f-score"], kind = "line", legend_loc = 4, directory = mtor_dir_results_graphs()):
    """ data - pandas data frame
        kind - bar or lines 
        Plots precision, recall, f-score, exports to file """
    if x_ticks is None:
        x_ticks = data["name"]
    import matplotlib
    import matplotlib.pyplot
    matplotlib.style.use('ggplot')
    
    matplotlib.pyplot.figure(); 
    data[columns].plot( kind = kind);
    matplotlib.pyplot.legend(loc = legend_loc)
    matplotlib.pyplot.xticks( range(len(data)), x_ticks, rotation = 30, ha ="right")
    matplotlib.pyplot.ylim((0,100))
    matplotlib.pyplot.savefig( "%s%s-%s.pdf" % (directory, name, kind))

def plot_comparison( datasets, dataset_names, name, x_ticks,
                             columns = ["precision", "recall", "f-score"], legend_loc = 4,
                             directory = mtor_dir_results_graphs()):
    """ Plot data side by side"""
    import matplotlib
    import matplotlib.pyplot
    matplotlib.style.use('ggplot')
    
    for column in columns:
        matplotlib.pyplot.figure();
        for d in datasets:
            matplotlib.pyplot.plot( d[column]);
        matplotlib.pyplot.legend( dataset_names, loc = legend_loc)
        matplotlib.pyplot.xticks( range(len(datasets[0])), x_ticks, rotation = 30, ha ="right")
        matplotlib.pyplot.ylim((0,100))
        matplotlib.pyplot.title( "%s %s" % (name, column))
        matplotlib.pyplot.savefig( "%s%s-%s-%s.pdf" % (directory,name, column, "-".join(dataset_names)))

########################################################################
########################################################################
# INITIALIZE

def initialize_mTORpathway_target():
    """ Initialize target.networkx.pickle """
    networkx_analysis.load_pathway( "TARGET", 
                                   input_file = mtor_dir() + "mTORpathway-sbml.xml", 
                                   output_file = mtor_dir() + "TARGET.networkx.pickle")

def initialize_mTORpathway_source( source = "DEFAULT"):
    """ initialize source NLP, ANN or others, None is default """
    # initialize annotation
    prefix = mtor_dir() + "events-" + source + "/"
    suffix = ".sbml.xml"
    files_target = glob.glob(  prefix + "*" + suffix)
    for f in files_target:
        print( "Processing %s" % f)
        id = f[len(prefix):-len(suffix)]
        name = "%s-%s" % (source,id)
        networkx_analysis.load_pathway( name, 
                                       input_file = f, 
                                       output_file = f + ".networkx.pickle",
                                       prefix = id + "_")
    print( "Combining graphs")
    files = glob.glob( prefix + "*.networkx.pickle")
    print( "Loading all graphs ...")
    graphs = [pickle.load( open( f, "rb")) for f in files]
    print( "Composing a single graph ...")
    graph = networkx.Graph()
    for g in graphs:
        graph = networkx.compose( graph, g, name = source)
    graph.name = source
    output_file = "%s/%s.networkx.pickle" % (mtor_dir(), source)
    print( "Exporting single graph to %s" % output_file)
    pickle.dump( graph, open( output_file, "wb"))


########################################################################
########################################################################
# run analysis

def run_simple_stats( dataset = "TARGET"):
    """ Export simple stats nodes, reaction numbers etc. """
    
    f = "%s%s.networkx.pickle" % (mtor_dir(), dataset)
    print( "%s: Processing %s from %s" % (now(), dataset,f))
    graph = networkx.read_gpickle( f)
    print( "%s:%s:%s: Filter isolated nodes" % (now(),dataset,graph))
    graph_no_isolated_nodes = networkx_analysis.filter_graph_remove_isolated_nodes( graph)
    
    for graph in [ graph, graph_no_isolated_nodes]:
        export_file = "%s%s-simple-stats.pickle" % (mtor_dir_results_statistics(), graph.name)
        networkx_analysis.run_analysis( graph, export_file)

def run_simple_stats_overlap( dataset_1 = "TARGET",
                              dataset_2 = "DEFAULT"):
                                  
    print( "%s: Processing %s/%s" % (now(), dataset_1, dataset_2))
    f1 = "%s%s.networkx.pickle" % (mtor_dir(), dataset_1)
    f2 = "%s%s.networkx.pickle" % (mtor_dir(), dataset_2)
    print( "Loading %s/%s" % (dataset_1, f1))
    graph_1 = networkx.read_gpickle( f1)
    print( "%s: Loading successful %s/%s/%s" % (now(), dataset_1, graph_1.name, f1))
    print( "%s: Loading %s/%s" % (now(), dataset_2, f2))
    graph_2 = networkx.read_gpickle( f2)
    print( "%s: Loading successful %s/%s/%s" % (now(), dataset_2, graph_2.name, f2))

    print( "%s: Computing %s-NO-ISOLATED-NODES" % (now(), dataset_2))
    graph_2_no_isolated_nodes = networkx_analysis.filter_graph_remove_isolated_nodes( graph_2)
    print( "%s: Finished computing %s-NO-ISOLATED-NODES/%s" % (now(), dataset_2, graph_2_no_isolated_nodes.name))
            
    ##################### SPECIES, REACTIONS, COMPARTMENTS OVERLAP
    networkx_analysis.run_analysis_signatures( graph_1, graph_2,
                                              export_file = "%s%s__%s-simple-stats-overlap.pickle" % (mtor_dir_results_statistics(), graph_1.name, graph_2.name))
    networkx_analysis.run_analysis_signatures( graph_1, graph_2_no_isolated_nodes,
                                              export_file = "%s%s__%s-simple-stats-overlap.pickle" % (mtor_dir_results_statistics(), graph_1.name, graph_2_no_isolated_nodes.name))
    

def run_node_overlap( node_match_fn_name,
                     dataset_1 = "TARGET",
                     dataset_2 = "ANN",
                     export_results = True,
                     export_results_prefix = mtor_dir_results() + "results-nodes-overlap-max",
                     compute_overlap_for_no_isolated_nodes = False,
                     ignore_existing = True):
    print( "Processing %s/%s %s" % ( dataset_1, dataset_2, node_match_fn_name))
    f1 = "%s%s.networkx.pickle" % (mtor_dir(), dataset_1)
    print( "Loading %s at %s" % (dataset_1, f1))
    graph_1 = networkx.read_gpickle( f1)
    f2 = "%s%s.networkx.pickle" % (mtor_dir(), dataset_2)
    print( "Loading %s at %s" % (dataset_2, f2))
    graph_2 = networkx.read_gpickle( f2)
    node_match_species = filter( lambda n: n.__name__ == node_match_fn_name, nodes_overlap_species_nm_fns)
    node_match_reaction = filter( lambda n: n.__name__ == node_match_fn_name, nodes_overlap_reactions_nm_fns)

    if compute_overlap_for_no_isolated_nodes:
        print( "Computing %s-NO-ISOLATED-NODES" % graph_2.name)
        graph_2_no_isolated_nodes = networkx_analysis.filter_graph_remove_isolated_nodes( graph_2)
    
    if node_match_species:
        node_match_species = node_match_species[0]
        networkx_analysis.run_analysis_nodes_overlap_max( networkx_analysis.filter_species( graph_1),
                                                         networkx_analysis.filter_species( graph_2),
                                                         node_match = node_match_species,
                                                         export_results = export_results,
                                                         export_results_prefix = export_results_prefix,
                                                         ignore_existing = ignore_existing);
        if compute_overlap_for_no_isolated_nodes:
            networkx_analysis.run_analysis_nodes_overlap_max( networkx_analysis.filter_species( graph_1),
                                                             networkx_analysis.filter_species( graph_2_no_isolated_nodes),
                                                             node_match = node_match_species,
                                                             export_results = export_results,
                                                             export_results_prefix = export_results_prefix,
                                                             ignore_existing = ignore_existing);
            
    if node_match_reaction:
        node_match_reaction = node_match_reaction[0]
        networkx_analysis.run_analysis_nodes_overlap_max( networkx_analysis.filter_reactions( graph_1),
                                                         networkx_analysis.filter_reactions( graph_2),
                                                         node_match = node_match_reaction,
                                                         export_results = export_results,
                                                         export_results_prefix = export_results_prefix,
                                                         ignore_existing = ignore_existing);
        if compute_overlap_for_no_isolated_nodes:
            networkx_analysis.run_analysis_nodes_overlap_max( networkx_analysis.filter_reactions( graph_1),
                                                             networkx_analysis.filter_reactions( graph_2_no_isolated_nodes),
                                                             node_match = node_match_reaction,
                                                             export_results = export_results,
                                                             export_results_prefix = export_results_prefix,
                                                             ignore_existing = ignore_existing);
    print( "Finished processing %s/%s %s" % ( dataset_1, dataset_2, node_match_fn_name))

def run_subgraphs_overlap( node_match_fn_name,
                          dataset_1 = "TARGET",
                          dataset_2 = "ANN",
                          export_results = True,
                          compute_overlap_for_with_isolated_nodes = False,
                          ignore_existing = True):
    print( "Processing %s/%s" % ( dataset_1, dataset_2))
    f1 = "%s%s.networkx.pickle" % (mtor_dir(), dataset_1)
    print( "Loading %s at %s" % (dataset_1, f1))
    graph_1 = networkx.read_gpickle( f1)
    f2 = "%s%s.networkx.pickle" % (mtor_dir(), dataset_2)
    print( "Loading %s at %s" % (dataset_2, f2))
    graph_2 = networkx.read_gpickle( f2)
    print( "Computing %s-NO-ISOLATED-NODES" % graph_2.name)
    graph_2_no_isolated_nodes = networkx_analysis.filter_graph_remove_isolated_nodes( graph_2)
    
    node_match = filter( lambda n: n.__name__ == node_match_fn_name, subgraphs_overlap_node_match_fns)[0]
    #### SUBGRAPH OVERLAP
    if compute_overlap_for_with_isolated_nodes:
        networkx_analysis.run_analysis_subgraphs_overlap( graph_1, graph_2,
                                       node_match = node_match, 
                                        export_results = export_results,
                                        export_results_prefix = mtor_dir_results() + "results-subgraphs-overlap-max",
                                        ignore_existing = ignore_existing)
                                            
    networkx_analysis.run_analysis_subgraphs_overlap( graph_1, graph_2_no_isolated_nodes,
                                   node_match = node_match, 
                                   export_results = export_results,
                                   export_results_prefix = mtor_dir_results() + "results-subgraphs-overlap-max",
                                   ignore_existing = ignore_existing)

    print( "Finished processing %s/%s %s" % ( dataset_1, dataset_2, node_match_fn_name))
                                                      
########################################################################
########################################################################
# node match plotting and printing

def print_AND_export_node_match_results( graph_1_name = "TARGET-SPECIES", 
                                        graph_2_name = "ANN-SPECIES", 
                                        node_match_names = ["nm_name_clean_approx", "nm_gene_id_intersect", "nm_bqbiol_is_equal", "nm_bqbiol_is_overlaps"]):
    import pickle
    for node_match_name in node_match_names:
        file_name = "results-nodes-overlap-max__%s__%s__%s.pickle" % (graph_1_name, graph_2_name, node_match_name)
        [graph_1, graph_2, matches] = pickle.load( open(file_name, "rb"))
        networkx_analysis.print_node_match_result( graph_1, graph_2, matches, node_match_name = node_match_name,
                                                  export_matches = file_name + ".txt")

# print_AND_export_node_match_results()
# print_AND_export_node_match_results( graph_1_name = "TARGET-SPECIES", graph_2_name = "NLP-SPECIES", node_match_names = ["nm_bqbiol_is_equal", "nm_bqbiol_is_overlaps"])

def node_match_results_to_pandas_dataframe( sources = ["ANN", "DEFAULT", "DEFAULT+GE11+MTOR"]):
    """ Exports species and reactions node match results as pandas dataframe """
    print("Loading reaction node match data")
    for source in sources:
        node_match_reaction_data = []
        for nm in nodes_overlap_reactions_nm_fns:
            graph_1_name = "TARGET-REACTIONS"
            graph_2_name = "%s-REACTIONS" % source
            f = "%s__%s__%s__%s.pickle" % (mtor_dir_results() + "results-nodes-overlap-max", graph_1_name, graph_2_name, nm.__name__)
            print( "Loading %s" % f)
            [graph_1, graph_2, matches] = pickle.load( open( f, "rb"))
            precision, recall, f_score = networkx_analysis.get_nodes_overlap_max_result_precision_recall_f_score( graph_1, graph_2, matches)
            node_match_reaction_data.append( {"target" : graph_1_name, "source" : graph_2_name, "name" : nm.__name__ , "precision" : precision, "recall" : recall, "f-score" : f_score})
        node_match_reaction_data = pandas.DataFrame(node_match_reaction_data)
        export_file_name = "%s%s-%s" % (mtor_dir_results_statistics(), graph_2_name, "node-match-reactions-result-statistics.pickle")
        print( "Exporting %s" % export_file_name)
        pandas.to_pickle( node_match_reaction_data, export_file_name)
        
    print("Loading species node match data")
    for source in sources:
        node_match_species_data = []
        for nm in nodes_overlap_species_nm_fns:
            graph_1_name = "TARGET-SPECIES"
            graph_2_name = "%s-SPECIES" % source
            f = "%s__%s__%s__%s.pickle" % (mtor_dir_results() + "results-nodes-overlap-max", graph_1_name, graph_2_name, nm.__name__)
            print( "Loading %s" % f)
            [graph_1, graph_2, matches] = pickle.load( open( f, "rb"))
            precision, recall, f_score = networkx_analysis.get_nodes_overlap_max_result_precision_recall_f_score( graph_1, graph_2, matches)
            node_match_species_data.append( {"target" : graph_1_name, "source" : graph_2_name, "name" : nm.__name__ , "precision" : precision, "recall" : recall, "f-score" : f_score})
        node_match_species_data = pandas.DataFrame(node_match_species_data)
        node_match_species_data[node_match_species_data["source"] == "ANN-SPECIES"]
        export_file_name = "%s%s-%s" % (mtor_dir_results_statistics(), graph_2_name, "node-match-species-result-statistics.pickle")
        print( "Exporting %s" % export_file_name)
        pandas.to_pickle( node_match_species_data, export_file_name)

#node_match_results_to_pandas_dataframe(sources = ["ANN", "DEFAULT", "DEFAULT+GE11+MTOR"])

def plot_node_match_results( source = "ANN"):
    """ exports line and bar plots of the statistics for node matching, requires pandas file of statistics """
    #node_match_results_to_pandas_dataframe()
    
    nm_species = pandas.read_pickle( mtor_dir_results() + "node_match_species_result_statistics.pickle")
    # filter data 
    filtered = numpy.logical_or.reduce( [nm_species["name"] == n.__name__ for n in nodes_overlap_species_nm_fns_paper])
    nm_species = nm_species[filtered].reset_index()
    nm_reactions = pandas.read_pickle( mtor_dir_results() + "node_match_reactions_result_statistics.pickle")
    
    plot_precision_recall_f_score(nm_species[nm_species["source"] == "%s-SPECIES" % source].reset_index(), "node-match-species-%s" % source, x_ticks = nodes_overlap_species_nm_fns_paper_names)
    plot_precision_recall_f_score(nm_reactions[nm_reactions["source"] == "%s-REACTIONS"].reset_index(), "node-match-reactions-%s" % source, x_ticks = nodes_overlap_reactions_nm_fns_names)
    plot_precision_recall_f_score(nm_species[nm_species["source"] == "%s-SPECIES" % source].reset_index(), "node-match-species-%s" % source, kind = "bar", x_ticks = nodes_overlap_species_nm_fns_paper_names)
    plot_precision_recall_f_score(nm_reactions[nm_reactions["source"] == "%s-REACTIONS"].reset_index(), "node-match-reactions-%s" % source, kind = "bar", x_ticks = nodes_overlap_reactions_nm_fns_names)

#plot_node_match_results("ANN")
#plot_node_match_results("DEFAULT")

def plot_node_match_results_comparison( sources = ["ANN","DEFAULT"]):
    """ plot node match comparisons for different sources """
    nm_species = pandas.read_pickle( mtor_dir_results() + "node_match_species_result_statistics.pickle")
    # filter data 
    filtered = numpy.logical_or.reduce( [nm_species["name"] == n.__name__ for n in nodes_overlap_species_nm_fns_paper])
    nm_species = nm_species[filtered].reset_index()
    nm_reactions = pandas.read_pickle( mtor_dir_results() + "node_match_reactions_result_statistics.pickle")

    datasets_species = [nm_species[nm_species["source"] == "%s-SPECIES" % source].reset_index() for source in sources]
    datasets_reactions = [nm_reactions[nm_reactions["source"] == "%s-REACTIONS" % source].reset_index() for source in sources]
    
    plot_comparison( datasets_species, sources, name = "node-match-species-comparison", x_ticks = nodes_overlap_species_nm_fns_paper_names)
    plot_comparison( datasets_reactions, sources, name = "node-match-reactions-comparison", x_ticks = nodes_overlap_reactions_nm_fns_names)
    
#plot_node_match_results_comparison()
    
########################################################################
########################################################################
# subgraph match plotting and printing

def subgraph_match_result_to_pandas_dataframe( target = "TARGET", source = "ANN"):
    """ Load subgraph data into pandas """
    print("%s:%s/%s: Load subgraph matching data" % (now(), target, source))
    
    edge_match = networkx_analysis.edge_match_exact
    results_prefix = mtor_dir_results() + "results-subgraphs-overlap-max"
    data = []
    for node_match in subgraphs_overlap_node_match_fns:
        f = "%s__%s__%s__%s__%s.pickle" % (results_prefix, target, source, node_match.__name__, edge_match.__name__)
        print("Loading %s/%s %s %s" % (target, source, node_match.__name__, f))
        [graph_1, graph_2, results_subgraphs, results_matches] = pickle.load( open( f, "rb"))
        
        result = networkx_analysis.get_subgraphs_overlap_max_results_precision_recall_f_score( graph_1 = graph_1, graph_2 = graph_2, results_subgraphs = results_subgraphs, results_matches = results_matches)
        result["target"] = target
        result["source"] = source
        result["node_match"] = node_match.__name__ 
        result["edge_match"] = edge_match.__name__ 
        data.append( result)

    data = pandas.DataFrame(data)
    export_file_name = "%ssubgraphs-overlap-max-statistics-%s__%s.pickle" % (mtor_dir_results_statistics(), target, source)
    print("%s:%s/%s: Exporting %s" % (now(), target, source, export_file_name))
    pandas.to_pickle( data, export_file_name)
    print("%s:%s/%s: Finished subgraph_match_result_to_pandas_dataframe" % (now(), target, source))

def subgraph_match_results_to_pandas_dataframe( sources = ["ANN", "DEFAULT", "DEFAULT+GE11+MTOR"]):
    """ Load subgraph data into pandas """
    print("Load subgraph matching data")
    target_name = "TARGET"
    edge_match = networkx_analysis.edge_match_exact
    results_prefix = mtor_dir_results() + "results-subgraphs-overlap-max"
    data = []
    for source_name in sources:
        for node_match in subgraphs_overlap_node_match_fns:
            f = "%s__%s__%s__%s__%s.pickle" % ( results_prefix, target_name, source_name, node_match.__name__, edge_match.__name__)
            print("Loading %s/%s %s %s" % (target_name,source_name, node_match.__name__, f))
            [graph_1, graph_2, results_subgraphs, results_matches] = pickle.load( open( f, "rb"))
            
            result = networkx_analysis.get_subgraphs_overlap_max_results_precision_recall_f_score( graph_1 = graph_1, graph_2 = graph_2, results_subgraphs = results_subgraphs, results_matches = results_matches)
            result["target"] = target_name
            result["source"] = source_name
            result["node_match"] = node_match.__name__ 
            result["edge_match"] = edge_match.__name__ 
            data.append( result)
    
    data = pandas.DataFrame(data)
    export_file_name = mtor_dir_results() + "results-subgraphs-overlap-max-statistics.pickle"
    print( "Exporting %s" % export_file_name)
    pandas.to_pickle( data, export_file_name)

# subgraph_match_results_to_pandas_dataframe( sources = ["ANN-NO-ISOLATED-NODES", "DEFAULT-NO-ISOLATED-NODES", "DEFAULT+GE11-NO-ISOLATED-NODES", "DEFAULT+GE11+MTOR-NO-ISOLATED-NODES", "DEFAULT+GE11+PC13+MTOR-NO-ISOLATED-NODES", "SK_DECISIONTREE+GE11+PC13+MTOR-NO-ISOLATED-NODES", "SK_MLP+GE11+MTOR-NO-ISOLATED-NODES", "SK_MULTINOMIALNB+GE11+MTOR-NO-ISOLATED-NODES", "SK_MULTINOMIALNB+GE11+PC13+MTOR-NO-ISOLATED-NODES"])

def plot_subgraph_results():
    """ Using pandas dataframes export all results
    Also prints some tables in latex format"""
    import collections, os, pandas, matplotlib, matplotlib.pyplot, tabulate
    matplotlib.style.use('ggplot')

    files = glob.glob( "%s/*-simple-stats.pickle" % mtor_dir_results_statistics())
    df = pandas.DataFrame([pandas.read_pickle(f) for f in files])
    # replace +MTOR with +ANN and DEFAULT with SVM
    for n in set(df["name"]):
        if "+MTOR" in n:
            df.replace( n, n.replace("+MTOR", "+ANN"), inplace = True)
            n = n.replace("+MTOR", "+ANN")
        if "DEFAULT" in n:
            df.replace( n, n.replace("DEFAULT", "SVM"), inplace = True)
            n = n.replace("DEFAULT", "SVM")
        for r1, r2 in [("SK_DECISIONTREE","DT"), ("SK_MULTINOMIALNB", "MNNB"), ("SK_MLP", "MLP"), ("SK_RF", "RF")]:
            if r1 in n:
                df.replace( n, n.replace( r1, r2), inplace = True)
                n = n.replace( r1, r2)
    
    # select all without -NO-ISOLATED-NODES
    _df = df[df["name"].apply( lambda x: not "-NO-ISOLATED-NODES" in x)]
    # remove ANN and DEFAULT
    _df = _df[_df["name"].apply( lambda x: not (x == "SVM" or x == "ANN"))]
    # replace SK_DECISIONTREE, DT, 
    for r1, r2 in [("SK_DECISIONTREE","DT"), ("SK_MULTINOMIALNB", "MNNB"), ("SK_MLP", "MLP"), ("SK_RF", "RF")]:
        _df.replace( r1, r2, inplace = True)
    _df = _df.reset_index()
    
    # individual graph statistics (nr nodes etc)
    columns = [["name", "# nodes", "# species", "# reactions"],
               ["name", "# nodes", "# species", "# reactions", "# compartments"],
               ["name", "# edges", "# edges reactant", "# edges product", "# edges modifier"],
               ["name", "isolates # nodes", "# isolated subgraphs"],
               ["name", '# isolated subgraphs', 'subgraphs # edges max', 'subgraphs # edges mean', 'subgraphs # edges median', 'subgraphs # edges min', 'subgraphs # nodes max', 'subgraphs # nodes mean', 'subgraphs # nodes median', 'subgraphs # nodes min', 'subgraphs # subgraphs'],
               ["name", "# species", "# reactions", "# compartments", "# edges", "# edges reactant", "# edges product", "# edges modifier", "isolates # nodes", "# isolated subgraphs"]]
    for c in columns:
        d = _df[c].set_index("name").sort_index()    
        matplotlib.pyplot.figure()
        d.plot(kind = "bar")
        matplotlib.pyplot.xticks( rotation = 30, horizontalalignment = "right")
        f = "%s%s.pdf" % (mtor_dir_results_graphs(), "-".join(c).replace(" ", "").replace("#","nr_"))
        try:
            f = f[:255]
        except:
            pass
        matplotlib.pyplot.savefig( f)
        print(tabulate.tabulate( d, headers=c, tablefmt = "latex"))

    input_files = glob.glob( "%s/subgraphs-overlap-max-statistics-*" % mtor_dir_results_statistics())
    graph_output_dir = mtor_dir_results_graphs()
    df = pandas.concat([pandas.read_pickle( f) for f in input_files])
    # replace node_match names with names in paper
    for nm in set(df["node_match"]):
        df.replace( nm, subgraphs_overlap_node_match_fns_names_map[nm], inplace = True)
    # select -NO-ISOLATED-NODES
    df = df[df["source"].apply( lambda x: "-NO-ISOLATED-NODES" in x)]
    # remove trailing NO-ISOLATED-NODES
    for s in set(df["source"]):
        if "-NO-ISOLATED-NODES" in s:
            df.replace( s, s.replace( "-NO-ISOLATED-NODES", ""), inplace = True)
            s = s.replace( "-NO-ISOLATED-NODES", "")
        if "DEFAULT" in s:
            df.replace( s, s.replace( "DEFAULT", "SVM"), inplace = True)
            s = s.replace( "DEFAULT", "SVM")
        for r1, r2 in [("SK_DECISIONTREE","DT"), ("SK_MULTINOMIALNB", "MNNB"), ("SK_MLP", "MLP"), ("SK_RF", "RF")]:
            if r1 in s:
                df.replace( s, s.replace( r1, r2), inplace = True)
                s = s.replace( r1, r2)
    # replace +MTOR with +ANN
    for n in set(df["source"]):
        if "+MTOR" in n:
            df.replace( n, n.replace("+MTOR", "+ANN"), inplace = True)
    # remove DEFAULT and ANN
    df = df[df["source"].apply( lambda x: not (x == "SVM" or x == "ANN"))]
    # replace SK_DECISIONTREE, DT, 
    df = df.reset_index()

    # plot complete summary
    for i in ["species", "reaction", "edge"]:
        for j in ["f-score", "precision", "recall"]:
            matplotlib.pyplot.figure()        
            
            df.pivot( "source", "node_match", "%s %s" % (i, j)).plot( kind = "bar")
            matplotlib.pyplot.ylim( (0, 100))
            matplotlib.pyplot.xticks( rotation = 30, horizontalalignment = "right")
            output_file = "%s/subgraphs-overlap-results-%s-%s.pdf" % (graph_output_dir, i, j)
            print( "Exporting %s " % output_file)        
            matplotlib.pyplot.savefig( output_file)

    # plot summary for each node-match
    columns = [ ["%s %s" % (i,j ) for j in ["precision", "recall", "f-score"]] for i in ["species", "reaction", "edge"]] + [["species f-score"],["reaction f-score"], ["edge f-score"]]
    columns += [ [ "%s %s" % (i,j) for i in ["species", "reaction", "edge"]] for j in ["precision", "recall", "f-score"]]
    for nm in set(df["node_match"]):
        d = df[df["node_match"] == nm].reset_index(drop=True).sort_index()
        for c in columns:
            matplotlib.pyplot.figure()
            d.pivot( "source", "node_match")[c].plot(kind = "bar")
            matplotlib.pyplot.ylim( (0, 100))        
            matplotlib.pyplot.xticks( rotation = 30, horizontalalignment = "right")
            output_file = "%ssubgraphs-overlap-results-%s---%s.pdf" % ( graph_output_dir, nm.replace("/","+").replace(", ","_"), "+".join(c).replace(" ", "-"))
            print( "Exporting %s " % output_file)        
            matplotlib.pyplot.savefig( output_file)
    
    # plot summary per data set
    for clf in set(df["source"]):
        d = df[df["source"] == clf].reset_index(drop=True).sort_index()
        for c in columns:
            matplotlib.pyplot.figure()
            d.pivot( "node_match", "source").reindex(index = subgraphs_overlap_node_match_fns_names)[c].plot(kind = "bar")
            matplotlib.pyplot.ylim( (0, 100))        
            matplotlib.pyplot.xticks( rotation = 30, horizontalalignment = "right")
            output_file = "%ssubgraphs-overlap-results-%s---%s.pdf" % ( graph_output_dir, clf, "+".join(c).replace(" ", "-"))
            print( "Exporting %s " % output_file)        
            matplotlib.pyplot.savefig( output_file)
            
    # boxplots for species/reactions/edge precision/recall/f-scpre
    for c in [ "%s %s" % (i,j ) for i in ["species", "reaction", "edge"] for j in ["precision", "recall", "f-score"]]:
        matplotlib.pyplot.figure()
        df[["source","node_match", c]].boxplot(by="source")
        matplotlib.pyplot.ylim( (0, 100))        
        matplotlib.pyplot.xticks( rotation = 30, horizontalalignment = "right")
        output_file = "%ssubgraphs-overlap-results-box-plot-%s.pdf" % ( graph_output_dir, c.replace( " ", "-"))
        print( "Exporting %s " % output_file)        
        matplotlib.pyplot.savefig( output_file)
    
        
    # tables for clf (f-score) for precision, recall and f-score separately
    for j in ["precision", "recall", "f-score"]:
        columns = [ "%s %s" % (i,j ) for i in ["species", "reaction", "edge"]]
        node_matches = set(df["node_match"])
        data = []
        for c in columns:
            for nm in node_matches:
                argmax = numpy.argmax( df[ df["node_match"] == nm][c])
                source = df.iloc[argmax]["source"]
                score = df.iloc[argmax][c]
                clf = source.split("+")[0]
                data.append( { "node_match": nm,
                               "column" : c,
                               "source" : source,
                               "score" : score,
                               "source_score" : "%s (%.0f)" % (source, score)})
        print( tabulate.tabulate( pandas.DataFrame( data).pivot( "node_match", "column", "source_score").reindex(index = subgraphs_overlap_node_match_fns_names), 
                                 headers = columns, tablefmt = "latex"))
    # meta-analysis of most often occuring best classifiers for precision, recall and f-score
    data = []
    for i in ["species", "reaction", "edge"]:
        for j in ["precision", "recall", "f-score"]:
            c = "%s %s" % (i,j)
            for nm in node_matches:
                argmax = numpy.argmax( df[ df["node_match"] == nm][c])
                source = df.iloc[argmax]["source"]
                score = df.iloc[argmax][c]
                data.append( {"node_match": nm,
                               "type1" : i,
                               "type2" : j,
                               "source" : source,
                               "score" : score})
    _df = pandas.DataFrame( data)
    for d in ["recall","precision","f-score"]:
        counter = collections.Counter(_df[_df["type2"] == d]["source"])
        matplotlib.pyplot.figure(); pandas.DataFrame( [{"source" : i, "count" : v} for i, v in counter.iteritems()]).set_index("source").sort_values( by="count")["count"].plot( kind = "bar")
        matplotlib.pyplot.title( "Best classifiers for %s histogram" % d)
        matplotlib.pyplot.xticks( rotation = 30, ha ="right")
        matplotlib.pyplot.ylim((0,60))
        output_file = "%ssubgraphs-overlap-results-best-classifiers-histogram-micro-%s.pdf" % ( graph_output_dir, d)
        print( "Exporting %s " % output_file)        
        matplotlib.pyplot.savefig( output_file)
        

    # compute macro recall/precision/f-score
    for p in ["precision", "recall"]:
        cols = [ "%s %s" % (s,p) for s in ["species","reaction","edge"]]
        df["macro %s" % p] = numpy.sum(df[cols], axis = 1) / float(len(cols))
    df["macro f-score"] = 2. * df["macro precision"] * df["macro recall"] / (df["macro precision"] + df["macro recall"])
    
    df = df[df["source"].apply( lambda x: not (x == "ANN"))]
    df_non_svm = df[df["source"].apply( lambda x: not (x == "SVM"))]
    
    this_year_max_fscore = [numpy.max( df_non_svm[df_non_svm["node_match"] == nm]["macro f-score"]) for nm in subgraphs_overlap_node_match_fns_names]
    last_year_f_score = df[df["source"] == "SVM"].set_index("node_match").reindex(subgraphs_overlap_node_match_fns_names)["macro f-score"]
    _df = pandas.DataFrame( { "max f-score" : this_year_max_fscore, "spranger2016 f-score" : last_year_f_score})
    print(tabulate.tabulate( _df, headers = _df.columns, tablefmt = "latex", floatfmt=".1f"))
    
    macro_columns = ["macro %s" % (j) for j in ["precision", "recall", "f-score"]]
    # boxplot macro recall/precision/f-score
    for c in macro_columns:
        matplotlib.pyplot.figure()
        df[["source","node_match", c]].boxplot(by="source")
        matplotlib.pyplot.ylim( (0, 100))
        matplotlib.pyplot.title( c)
        matplotlib.pyplot.xticks( rotation = 90, horizontalalignment = "right")
        output_file = "%ssubgraphs-overlap-results-box-plot-%s.pdf" % ( graph_output_dir, c.replace( " ", "-"))
        print( "Exporting %s " % output_file)        
        matplotlib.pyplot.savefig( output_file)
    
    # plot node_match vs source value: macro precision/recall/f-score
    # this is for all 
    for c in macro_columns:
        print("Table %s all results" % c)
        _df = df.pivot( "node_match", "source", c).reindex( index = subgraphs_overlap_node_match_fns_names)
        print(tabulate.tabulate( _df, tablefmt = "latex", headers = _df.columns, floatfmt=".1f"))
    
    # compute best classifiers per node_match and macro precision/recall/f-score
    data = []
    for nm in set(df["node_match"]):
        d = {"node_match": nm}
        for c in macro_columns:
            argmax = numpy.argmax( df[ df["node_match"] == nm][c])
            source = df.iloc[argmax]["source"]
            score = df.iloc[argmax][c]
            d["%s source" % c] = source
            d["%s score" % c] = score
            d[c] = "%s (%.f)" % (source,score)
        data.append( d)
    _df = pandas.DataFrame( data) # .set_index("node_match").reindex(index = subgraphs_overlap_node_match_fns_names)
    
    # plot histograms of best classifier
    for c in macro_columns:
        counter = collections.Counter(_df["%s source" % c])
        matplotlib.pyplot.figure(); 
        pandas.DataFrame( [{"source" : i, "count" : v} for i, v in counter.iteritems()]).set_index("source").sort_values( by="count")["count"].plot( kind = "bar")
        matplotlib.pyplot.title( "Best classifiers for macro %s histogram" % c)
        matplotlib.pyplot.xticks( rotation = 30, ha ="right")
        matplotlib.pyplot.ylim((0,20))
        output_file = "%ssubgraphs-overlap-results-best-classifiers-histogram-%s.pdf" % ( graph_output_dir, c.replace(" ", "-"))
        print( "Exporting %s " % output_file)        
        matplotlib.pyplot.savefig( output_file)
    
    ## print table node_match vs macro precision/recall/f-score, value is best source and score
    # compute best classifiers per node_match and macro precision/recall/f-score
    print( tabulate.tabulate( _df.set_index( "node_match").reindex(index = subgraphs_overlap_node_match_fns_names)[macro_columns], headers = columns, tablefmt = "latex"))
    figure = matplotlib.pyplot.figure(figsize=(10, 4))
    figure.add_subplot(111)
    _df.set_index( "node_match").reindex(index = subgraphs_overlap_node_match_fns_names)[["%s score" % m for m in macro_columns]].plot(kind="bar")
    matplotlib.pyplot.ylim((0,100))
    matplotlib.pyplot.xticks( rotation = 90, ha ="right")
    output_file = "%ssubgraphs-overlap-results-best-macro-precision-macro-recall-macro-f-score-per-node-match.pdf" % ( graph_output_dir)
    print( "Exporting %s " % output_file)        
    matplotlib.pyplot.savefig( output_file)

#plot_subgraph_results()

########################################################################
########################################################################
# MAIN

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--initialize-target',
                    action = "store_true",
                    dest = "initialize_target",
                    default = False,
                    help = "Initialize TARGET dataset")
    parser.add_argument( '--initialize-graph',
                    action = "store",
                    dest = "initialize_graph",
                    default = None,
                    help = "Initialize some dataset (e.g. DEFAULT, ANN etc)")
    parser.add_argument( '--simple-stats',
                    action = "store",
                    dest = "simple_stats",
                    default = None, # e.g. TARGET, DEFAULT, ANN
                    help = "Run a simple analysis on some dataset (e.g. TARGET, DEFAULT, ANN)")
    parser.add_argument( '--simple-stats-overlap',
                    action = "store",
                    dest = "simple_stats_overlap",
                    default = None, # e.g. DEFAULT, ANN
                    help = "Run a simple overlap analysis on TARGET vs parameter (e.g. DEFAULT, ANN)")
    parser.add_argument( '--ignore-existing',
                    action = "store_true",
                    dest = "ignore_existing",
                    default = False,
                    help = "Ignore existing results")
    ### NODE OVERLAP
    parser.add_argument( '--node-match-results-to-pandas-dataframe',
                    action = "store",
                    dest = "node_match_results_to_pandas_dataframe",
                    default = None,
                    help = "Translates networkx pickle results in to pandas data frame")
    parser.add_argument( '--print-node-overlap-fns',
                    action = "store_true",
                    dest = "print_node_overlap_fns",
                    default = False,
                    help = "Prints all current node overlap matching strategies")
    parser.add_argument( '--node-overlap-node-match-fn',
                    action = "store",
                    dest = "node_overlap_node_match_fn",
                    default = nodes_overlap_species_nm_fns[0].__name__, # e.g. 
                    help = "Node match function used in node overlap run")
    parser.add_argument( '--node-overlap',
                    action = "store",
                    dest = "node_overlap",
                    default = None,
                    help = "Run a node overlap analysis for TARGET vs [parameter] (e.g. DEFAULT, ANN)")
    ### SUBGRAPH OVERLAP
    parser.add_argument( '--print-subgraphs-overlap-fns',
                    action = "store_true",
                    dest = "print_subgraphs_overlap_fns",
                    default = False,
                    help = "Prints all current subgraph matching strategies")
    parser.add_argument( '--subgraphs-overlap-node-match-fn',
                    action = "store",
                    dest = "subgraphs_overlap_node_match_fn",
                    default = subgraphs_overlap_node_match_fns[0].__name__, # e.g. 
                    help = "Node match function used in subgraphs overlap run")
    parser.add_argument( '--subgraphs-overlap',
                    action = "store",
                    dest = "subgraphs_overlap",
                    default = None, # e.g. 
                    help = "Run a subgraph analysis for TARGET vs [parameter] (e.g. DEFAULT, ANN)")
    parser.add_argument( '--print-subgraphs-overlap-statistics',
                    action = "store",
                    dest = "print_subgraphs_overlap_statistics",
                    default = None, # e.g. 
                    help = "Print precision recall f-score (uses stored data) for TARGET vs [parameter] (e.g. DEFAULT, ANN)")
    parser.add_argument( '--subgraph-match-result-to-pandas-dataframe',
                    action = "store",
                    dest = "subgraph_match_result_to_pandas_dataframe",
                    default = None, # e.g. 
                    help = "Load all subgraph match results and create pandas dataframe")
    parser.add_argument( '--subgraphs-overlap-sbml',
                    action = "store",
                    dest = "subgraphs_overlap_sbml",
                    default = None, # e.g. 
                    help = "Export colored SBML xml files for all analyses (uses stored data) for TARGET vs [parameter] (e.g. DEFAULT, ANN)")
    parser.add_argument( '--helper-load-subgraph-result-export-to-pandas',
                    action = "store",
                    dest = "helper_load_subgraph_result_export_to_pandas",
                    default = None, # e.g. 
                    help = "Export panda files for subgraph result (pass node match fn name)")
    parser.add_argument( '--subgraphs-overlap-result-graphs',
                    action = "store_true",
                    dest = "subgraphs_overlap_result_graphs",
                    default = None, # e.g. 
                    help = "Export precision, recall, f-score (bar) graphs; also prints TEX table information")
    
    cmd = parser.parse_args()

    
    if cmd.initialize_target:
        networkx_analysis.initialize()
        initialize_mTORpathway_target()
    if cmd.initialize_graph:
        networkx_analysis.initialize()
        initialize_mTORpathway_source(cmd.initialize_graph)
    if cmd.simple_stats:
        networkx_analysis.initialize()
        run_simple_stats(cmd.simple_stats)
    if cmd.simple_stats_overlap:
        networkx_analysis.initialize()
        run_simple_stats_overlap( dataset_2 = cmd.simple_stats_overlap)
    
    if cmd.node_match_results_to_pandas_dataframe:
        node_match_results_to_pandas_dataframe(sources = [cmd.node_match_results_to_pandas_dataframe])
    if cmd.print_node_overlap_fns:
        print("\n".join( [fn.__name__ for fn in nodes_overlap_species_nm_fns]))
        print("\n".join( [fn.__name__ for fn in nodes_overlap_reactions_nm_fns]))
    if cmd.node_overlap:
        networkx_analysis.initialize()
        run_node_overlap( cmd.node_overlap_node_match_fn, 
                         dataset_2 = cmd.node_overlap,
                         ignore_existing = cmd.ignore_existing)

    if cmd.print_subgraphs_overlap_fns:
        print("\n".join( [fn.__name__ for fn in subgraphs_overlap_node_match_fns]))
    if cmd.subgraphs_overlap:
        networkx_analysis.initialize()
        run_subgraphs_overlap( cmd.subgraphs_overlap_node_match_fn, 
                              dataset_2 = cmd.subgraphs_overlap,
                              ignore_existing = cmd.ignore_existing)
    if cmd.print_subgraphs_overlap_statistics:
        for node_match in subgraphs_overlap_node_match_fns:
            networkx_analysis.print_analysis_subgraphs_overlap_results_from_file( "TARGET", 
                                                                                 cmd.print_subgraphs_overlap_statistics, 
                                                                                 node_match,
                                                                                 prefix = mtor_dir_results() + "results-subgraphs-overlap-max")
    if cmd.subgraph_match_result_to_pandas_dataframe:
        subgraph_match_result_to_pandas_dataframe( source = cmd.subgraph_match_result_to_pandas_dataframe)
    if cmd.subgraphs_overlap_sbml:
        edge_match = networkx_analysis.edge_match_exact
        prefix = mtor_dir() + "results/results-subgraphs-overlap-max"
        for isolated in ["","-NO-ISOLATED-NODES"]:
            graph_2_name = "%s%s" % (cmd.subgraphs_overlap_sbml,isolated)
            for node_match in subgraphs_overlap_node_match_fns:
                print( "%s Loading %s/%s" %  (now(), graph_2_name, node_match.__name__))
                [graph_1, graph_2, results_subgraphs, results_matches] \
                = pickle.load( open( "%s__%s__%s__%s__%s.pickle" % ( prefix, "TARGET", graph_2_name, node_match.__name__, edge_match.__name__), "rb"))
                print( "%s Computing results %s/%s" %  (now(), graph_2_name, node_match.__name__))
                matches_nodes_1, matches_nodes_2, \
                matches_edges_1, matches_edges_2, \
                species_1_matches, species_2_matches, \
                reactions_1_matches, reactions_2_matches = \
                networkx_analysis.get_subgraphs_overlap_max_results( graph_1, graph_2, results_subgraphs, results_matches)
                print( "%s Exporting SBML/XML %s/%s" %  (now(), graph_2_name, node_match.__name__))
                networkx_analysis.subgraph_overlaps_sbml( graph_1, matches_nodes_1, matches_edges_1,
                                        inn = mtor_dir() + 'mTORpathway-celldesigner.xml',
                                        out = mtor_dir() + 'results-sbml/mTORpathway-celldesigner-%s-%s-%s.xml' % (graph_2_name, node_match.__name__, edge_match.__name__),
                                        background_color_bounds = "00000000",
                                        background_color_lines = "90000000",
                                        matched_color = "FF00FF00")
    if cmd.helper_load_subgraph_result_export_to_pandas:
        print("Exporting subgraph matching data to pandas dataframe")
        print("Load subgraph matching data")
        target_name = "TARGET"
        source_name = "ANN"
        edge_match = networkx_analysis.edge_match_exact
        results_prefix = mtor_dir_results() + "results-subgraphs-overlap-max"
        node_match_name = cmd.helper_load_subgraph_result_export_to_pandas
        for fn in subgraphs_overlap_node_match_fns:
            node_match_name = fn.__name__
            for source_name in [ "%s%s" % (cmd.helper_load_subgraph_result_export_to_pandas, isolated) for isolated in ["", "-NO-ISOLATED-NODES"]]:
                print("Loading %s %s" % (source_name, node_match_name))
                [graph_1, graph_2, results_subgraphs, results_matches] = pickle.load( open( "%s__%s__%s__%s__%s.pickle" % ( results_prefix, target_name, source_name, node_match_name, edge_match.__name__), "rb"))
                
                result = networkx_analysis.get_subgraphs_overlap_max_results_precision_recall_f_score( graph_1 = graph_1, graph_2 = graph_2, results_subgraphs = results_subgraphs, results_matches = results_matches)
                result["target"] = target_name
                result["source"] = source_name
                result["node_match"] = node_match_name
                result["edge_match"] = edge_match.__name__ 
                export_file_name = mtor_dir_results() + "results-subgraphs-overlap-max-statistics-%s-%s.pickle" % (source_name, node_match_name)
                print( "Exporting %s" % export_file_name)
                pickle.dump( result, open( export_file_name, "wb"))
    if cmd.subgraphs_overlap_result_graphs:
        print( "Exporting subgraph overlap result (bar) graphs for data sets, also print TEX table information")
        f = mtor_dir_results() + "results-subgraphs-overlap-max-statistics.pickle"
        print( "Loading %s" % f)
        df = pandas.read_pickle( f)
        sources = set(df["source"])
        for source in sources:
            data = df[ df["source"] == source].sort_values( "node_match").reset_index()
            # sort by order in list of functions
            subgraphs_overlap_node_match_fns_real_names = [n.__name__ for n in subgraphs_overlap_node_match_fns]
            data["order"] = data["node_match"].apply( lambda v: subgraphs_overlap_node_match_fns_real_names.index( v))
            data = data.sort_values("order").reset_index()        
            for t in ["precision","recall","f-score"]:
                name = "subgraph-overlap-match-%s-%s-wo-constituents" % (source, t)
                columns = ["%s %s" % (s,t) for s in ["species", "reaction","edge"]]
                plot_precision_recall_f_score( data[:12], name, x_ticks = subgraphs_overlap_node_match_fns_names[:12], columns = columns, kind = "line", legend_loc = 2,
                                              directory = mtor_dir_results_graphs())
                plot_precision_recall_f_score( data[:12], name, x_ticks = subgraphs_overlap_node_match_fns_names[:12], columns = columns, kind = "bar", legend_loc = 2,
                                              directory = mtor_dir_results_graphs())
                name = "subgraph-overlap-match-%s-%s-w-constituents" % (source, t)
                plot_precision_recall_f_score( data[12:].reset_index(drop=True), name, x_ticks = subgraphs_overlap_node_match_fns_names[12:], columns = columns, kind = "line", legend_loc = 2,
                                              directory = mtor_dir_results_graphs())
                plot_precision_recall_f_score( data[12:].reset_index(drop=True), name, x_ticks = subgraphs_overlap_node_match_fns_names[12:], columns = columns, kind = "bar", legend_loc = 2,
                                              directory = mtor_dir_results_graphs())
        
            print( "----------- TEX BEGIN %s" % source)
            # export to TEX TABLE
            columns = [ "%s %s" % (i,j) for i in ["node","species","reaction","edge"] for j in ["precision","recall","f-score"]]
            for i, t in zip( range(len( data)), subgraphs_overlap_node_match_fns_names):
                s = t
                for column in columns:
                    s += " & %.2f" % data.iloc[i][column]
                s += "\\\\\\hline"
                print( s)
            print( "----------- TEX END %s" % source)
