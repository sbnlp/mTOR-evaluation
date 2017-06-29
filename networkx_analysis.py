import collections
import csv
import datetime
import fuzzywuzzy.fuzz
import fuzzywuzzy.process
import itertools
import joblib
import libsbml
import lxml
import lxml.etree
import networkx
import numpy
import os
import operator
import pickle
import re
import simstring
import sys

########################################################################
########################################################################
# Globals

# gene_map
GENE_MAP = None
# simstring
SIMSTRING_DB = None

SBO_NODES = None
#SBO_NODES = convert_xml_to_sbonodes()

########################################################################
########################################################################

def now():
    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

########################################################################
########################################################################

def exists( x, elements, test = lambda x,y : x == y):
    for y in elements:
        if test( x, y):
            return True
    return False

########################################################################
########################################################################
# remove_prefixes

PREFIXES = [ "acetylated ", "activated ", "associated ", \
    "bound ", \
    "catabolized ", "catalyzed ", "converted ", \
    "deacetylated ", "degradated ", "demethylated ", "dephosporylated ", "deubiquinated ", "dissociated ","deactivated ", \
    "expressed ", \
    "methylated ", \
    "positively ",\
    "negatively ", \
    "regulated ",\
    "phosphorylated ", 
    "regulated ",\
    "transcribed ", "translated ", \
    "ubiquitinated "]
    
def remove_prefixes( name):
    global PREFIXES
    new_name = name
    for prefix in PREFIXES:
        if prefix != None:
            new_name = new_name.replace( prefix, "")
    return new_name.strip()

########################################################################
########################################################################

def compute_all_is_a( node, nodes):
    all_parents = set( node["is_a"])
    for parent_id in node["is_a"]:
        all_parents.update( compute_all_is_a( nodes[parent_id], nodes))
    return all_parents
    
def convert_xml_to_sbonodes( file_name = "sbo.xml", output_file_name = "sbo.pickle"):

    # load nodes
    nodes = {}
    sbo_xml = lxml.etree.fromstring( open( file_name, "rt").read())
    
    for term in sbo_xml.xpath( "/*[local-name()='sbo']/*[local-name()='Term']"):
        id = term.find( "{http://www.biomodels.net/sbo}id").text
        name = term.find( "{http://www.biomodels.net/sbo}name").text
        is_a = [];
        if term.find( "{http://www.biomodels.net/sbo}is_a") is not None:
            is_a = [el.text for el in term.findall( "{http://www.biomodels.net/sbo}is_a")]
        nodes[id] = { "id" : id, "name" : name , "is_a" : is_a }
    
    # compute all is_a for fast lookup
    is_a_all = {}
    for node in nodes.itervalues():
        is_a_all[node["id"]] = compute_all_is_a( node, nodes)
    for node in nodes.itervalues():
        node["is_a"] = is_a_all[node["id"]]
    if output_file_name is not None:
        pickle.dump( nodes, open( output_file_name, "wb"))
    return nodes;


def sbo_is_a( sbo_1, sbo_2):
    "return true if sbo_1 is_a sbo_2 (if any of them is None, return true)"
    global SBO_NODES
    if sbo_1 == sbo_2 or sbo_1 == None or sbo_2 == None:
        return True
    elif sbo_1 in SBO_NODES:
        return sbo_2 in SBO_NODES[sbo_1]["is_a"];
    else:
        return False

def sbo_is_a2( sbo_1, sbo_2):
    "Return true if is a either direction"
    return sbo_is_a( sbo_1, sbo_2) or sbo_is_a( sbo_2, sbo_1)

def sbo_name( sbo_1):
    global SBO_NODES
    return SBO_NODES[sbo_1]["name"]

def load_sbo( file_name = "sbo.pickle"): 
    global SBO_NODES
    SBO_NODES = pickle.load( open( file_name, "rb"))

def sbo_export_graph():
    global SBO_NODES
    sbo_graph = networkx.DiGraph()
    for node in SBO_NODES:
        sbo_graph.add_node( node)
    for node in SBO_NODES.values():
        for parent in node["is_a"]:
            sbo_graph.add_edge( node["id"], parent)
    export_all_graph( sbo_graph, "sbo_graph")

def sbo_export_graph_nodes( nodes, file_prefix = "test"):
    """ exports hierarchy for SBO nodes"""
    global SBO_NODES
    sbo_graph = networkx.DiGraph()
    
    all_nodes = nodes + [ parent for n in nodes for parent in compute_all_is_a( n) ]
    for node in all_nodes:
        sbo_graph.add_node( node)
    
    for node in all_nodes:
        for parent in node["is_a"]:
            sbo_graph.add_edge( node["id"], parent)
    export_all_graph( sbo_graph, file_prefix)

def get_terms( miriam_urns):
    """ takes a list of miriam encoded urn, e.g. ['urn:miriam:GO:0016579', 'urn:miriam:SBO:0000330']
        and returns the strings ["GO:0016579", "SBO:0000330"] """
    return [ i[11:]for i in miriam_urns]

def get_sbo_terms( miriam_urns):
    """ takes a list of miriam encoded urn, e.g. ['urn:miriam:GO:0016579', 'urn:miriam:SBO:0000330']
        and returns the strings ["SBO:0000330"] """
    return [ i[11:]for i in miriam_urns if i.startswith( "urn:miriam:SBO:")]

def get_sbo_int( miriam_urns):
    """ takes a list of miriam encoded urn, e.g. ['urn:miriam:GO:0016579', 'urn:miriam:SBO:0000330']
        and returns the integers [330] """
    return [ int( i[15:]) for i in miriam_urns if i.startswith( "urn:miriam:SBO:")]

########################################################################
########################################################################

ST_SBO_GO_MAP = {  # degradation
    'acetylation': 'SBO:0000215',
    'activation': 'SBO:0000170',
    'association': 'SBO:0000297',
    'binding': 'SBO:0000297',
    'catabolism': 'GO:0009056',
    'catalysis': 'SBO:0000172',
    'conversion': 'SBO:0000182',
    'deacetylation': 'GO:0006476',
    'degradation': 'SBO:0000179',
    'demethylation': 'GO:0006482',
    'dephosphorylation': 'SBO:0000330',
    'deubiquitination': 'GO:0016579',
    'dissociation': 'SBO:0000180',
    'gene_expression': 'SBO:0000205',
    'inactivation': 'SBO:0000169',
    'localization': 'GO:0051179',
    'methylation': 'SBO:0000214',
    'negative_regulation': 'SBO:0000169',
    'pathway': 'SBO:0000375',
    'phosphorylation': 'SBO:0000216',
    'positive_regulation': 'SBO:0000170',
    'protein_catabolism': 'SBO:0000179',
    'regulation': 'SBO:0000168',
    'transcription': 'SBO:0000183',
    'translation': 'SBO:0000184',
    'transport': 'SBO:0000185',
    'ubiquitination': 'SBO:0000224'}

SBO_GO_ST_MAP = { v : k for k, v in ST_SBO_GO_MAP.iteritems()}

def sbo_go_name( urn_miriam):
    if urn_miriam.startswith( "urn:miriam:"):
        urn_miriam = urn_miriam[11:]
    if urn_miriam in SBO_GO_ST_MAP:
        return SBO_GO_ST_MAP[urn_miriam]
    elif urn_miriam.startswith( "SBO:"):
        return sbo_name( urn_miriam)
    else:
        return urn_miriam

def sbo_go_name_known( urn_miriam):
    if urn_miriam.startswith( "urn:miriam:"):
        urn_miriam = urn_miriam[11:]
    if urn_miriam in SBO_GO_ST_MAP:
        return True
    elif urn_miriam.startswith( "SBO:"):
        return True
    else:
        return False
########################################################################
########################################################################

def clean_name( name):
    return remove_prefixes( name.lower()).strip()

def clean_name2( name):
    return re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( name.lower())).strip()

def names( graph):
    return [graph.node[n].get("name") for n in graph.nodes() if graph.node[n].get("name")]

def names_clean( graph):
    return [ remove_prefixes( graph.node[n].get("name").lower()) for n in graph.nodes() if graph.node[n].get("name")]

def names_clean2( graph):
    return [ re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( graph.node[n].get("name").lower())) for n in graph.nodes() if graph.node[n].get("name")]

########################################################################
########################################################################

def sort_edge_signature( signature, graph):
    """ takes (species122,reaction122,"product") and returns (reaction122,species122,"product") """
    if signature[2] == "reactant" and graph.node[signature[0]]["type"] != "species":
        return (signature[1],signature[0],signature[2])
    elif signature[2] == "product" and graph.node[signature[1]]["type"] != "species":
        return (signature[1],signature[0],signature[2])
    elif signature[2] == "modifier" and graph.node[signature[0]]["type"] != "species":
        return (signature[1],signature[0],signature[2])
    else:
        return signature

def edge_signatures( graph):
    signatures = set([ sort_edge_signature( (edge[0], edge[1], edge[2]["type"]), graph) for edge in graph.edges( data = True)])
    assert( len(signatures) == len( graph.edges()))
    return signatures

########################################################################
########################################################################

def create_gene_map( chilibot = True, hugo = True, human_entrez = False):
    lists = []
    print( "create_gene_map")   
    print("Loading data")
    if chilibot:
        with open( "gene_list_chilibot.txt", "rt") as f:
            txt = f.read()
            for line in txt.strip().split("\n"):
                line = line.strip(";")
                synonyms = [ line.split( "|")[0].strip()] + line.split( "|")[1].split( ";")
                lists.append( set( [s.lower() for s in synonyms]))
    if hugo:
        with open('gene_list_hugo.txt', 'rU') as f:
            csv_list = csv.reader( f, delimiter = '\t')
            for row in csv_list:
                lists.append( set( [ s.lower() for s in filter( bool, row) if s != ""]))
    if human_entrez:
        with open('gene_list_human_entrez.txt', 'r') as f:
            lines = f.read().split("\n")
            lines.pop(0) # remove first line
            for line in lines:
                synonyms = [s.lower() for s in line.strip().split("\t")]
                synonyms.pop(0)
                lists.append( set(synonyms))
    
    print("Merging lists")
    
    dict_forward = {} # maps el : value
    dict_backward = {} # maps val : list of elements
    new_value_counter = 0
    
    for idx, l in enumerate(lists):
        if idx % 10000 == 0:
            print( "Processed %i" % idx)
        new_value_counter += 1 
        new_value = new_value_counter
        # compute overlap_values - those values overlapping
        overlap_values = set()
        for e in l:
            if e in dict_forward:
                overlap_values.add( dict_forward[e])
        elements = set(l) # initialize elements with known values
        if overlap_values != set():
            new_value = new_value_counter
            new_value_counter += 1 
            # update elements with known values
            for val in overlap_values:
                elements.update( dict_backward[val])
            # update dict_forward
            for e in elements:
                dict_forward[e] = new_value
            # update dict_backward
            for val in overlap_values:
                del dict_backward[val]
            dict_backward[new_value] = elements
        else: # no overlap found, just add elements to dicts
            for e in elements:
                dict_forward[e] = new_value
            dict_backward[new_value] = elements
                
    lists = list(dict_backward.values())            
    print("Merging lists finished (%i total sets)" % len( lists))
    
    print("Computing gene map")
    gene_map = {}
    for l in lists:
        listt = [ re.sub('[^a-zA-Z0-9-]', ' ', e.lower()) for e in l if e != ""]
        if listt != []:
            val = listt[0]
            for l in listt:
                gene_map[l] = val
    print("Computing gene map (%i total names/genes)" % len( gene_map))
    
    print("Exporting gene map")
    pickle.dump( gene_map, open( "gene_map.pickle", "wb"))
    
    return gene_map

def create_simstring_txt( gene_map):
    """ Creates gene_list.txt for usage in simstring db
        use: simstring -b -d gene_list.simstring < gene_list.txt
        afterwards to create simstring"""
    print( "create_simstring_txt")   
    with open( "gene_list.txt", "wt") as f:
        f.write( "\n".join( gene_map.keys() + list( set( gene_map.values()))))

def create_simstring_db():
    """ Creates simstring database
        use: simstring -b -d gene_list.simstring < gene_list.txt"""
    import commands
    print( "create_simstring_db")
    ret = commands.getstatusoutput('simstring -b -d gene_list.simstring < gene_list.txt')
    print( ret)
    print( "create_simstring_db finished")

def create_gene_map_AND_simstring_db():
    gene_map = create_gene_map()
    # gene_map = pickle.load( open( "gene_map.pickle", "rb"))
    create_simstring_txt( gene_map)
    create_simstring_db()

#######################

def map_gene_fuzzywuzzy( name, threshold = 90):
    global GENE_MAP
    assert(GENE_MAP)
    clean_name = clean_name2( name)
    if GENE_MAP.get( clean_name):
        return set( [GENE_MAP[clean_name]])
    else:
        results = set()
        for k in GENE_MAP.keys():
            if fuzzywuzzy.fuzz.ratio( clean_name, k) > threshold:
                results.add( GENE_MAP[k])
        if results != set():
            return results
        else:
            return None

def map_gene_simstring( name):
    "retrieves gene_map results by simstring matching and lookup"
    global GENE_MAP, SIMSTRING_DB
    assert( GENE_MAP and SIMSTRING_DB)
    clean_name = clean_name2( name)
    if GENE_MAP.get( clean_name):
        return set( [GENE_MAP[clean_name]])
    else:
        results = SIMSTRING_DB.retrieve( clean_name)
        if results:
            return set( [GENE_MAP[r] for r in results])
        else:
            return None

def export_mapping( mapping, file_name):
    with open( file_name, "wt") as f:
        f.write( "\n".join( [ "{} : {}".format( k, ",".join( [str(v) for v in values])) for k, values in mapping.itervalues()]))
    
def compute_simstring_coverage( names, thresholds = [ i/10.0 for i in range(1, 10)], measure = simstring.cosine):
    results = []
    for t in thresholds:
        db = simstring.reader( 'gene_list.simstring')
        db.measure = measure
        db.threshold = t
        results.append( [ True for n in names if map_gene_simstring(n, db)].count( True) / float( len( names)))
    return results
    
########################################################################
########################################################################

def export_graph( graph, graph_name, prog = "dot"):
    agraph = networkx.to_agraph( graph)
    ## "neato"|"dot"|"twopi"|"circo"|"fdp"|"nop"
    agraph.layout( prog = prog)
    file_name = graph_name + "_" + prog + ".pdf"
    agraph.draw( file_name)
    print( "Exported {}".format( file_name))

def export_all_graph( graph, graph_name):
    for prog in ["neato", "dot", "twopi", "circo", "fdp"]:
        export_graph( graph, graph_name, prog = prog)

########################################################################
########################################################################

def load_sbml( file_name):
    reader = libsbml.SBMLReader()
    document = reader.readSBML( file_name)
    print( "Loaded {} ({} errors)".format( file_name, document.getNumErrors()))
    return document

def get_participants_species( species, prefix, model):
    """ Takes an SBML species and returns its participants (mTOR)"""
    annotation = species.getAnnotation()
    if annotation == None:
        return []
    # retrieve path
    annotation_path_names = [ 'RDF', 'Participants']
    current_state = annotation
    for name in annotation_path_names:
        last_state = current_state
        current_state = None
        for i in xrange( last_state.getNumChildren()):
            if last_state.getChild(i).getName() == name:
                current_state = last_state.getChild(i)
        if current_state == None:
            break
    # retrieve participants
    participants = []
    if current_state != None:
        for idx in range( current_state.getNumChildren()):
            child = current_state.getChild( idx)
            if child.getName() != 'Participant':
                sys.stderr.write( "\nERROR: unexpected participant xml name {}".format( prefix + species.getId()))
                sys.stderr.flush()
            elif child.getAttrValue("participant") == "":
                sys.stderr.write( "\nERROR: unexpected participant attribute value {}".format( prefix + species.getId()))
                sys.stderr.flush()
            elif model.getSpecies( child.getAttrValue("participant")) == None:
                sys.stderr.write( "\nERROR: participant {} does not exist in model (species: {})".format( child.getAttrValue("participant"), prefix + species.getId()))
                sys.stderr.flush()
            else:
                participants.append( child.getAttrValue("participant"))
    return participants

def create_graph( model, prefix = "", ignore_participant_graph = False,
                 skip_uris = ["urn:miriam:reactome", "urn:miriam:pubmed", "urn:miriam:ec"]):

    graph = networkx.Graph();
    
    # add species
    for species in model.getListOfSpecies():
        bqbiol_is = []
        bqbiol_has_part = []
        bqbiol_has_version = []
        if species.getCVTerms() != None:
            for term in species.getCVTerms():
                uris = [ term.getResourceURI( idx) for idx in xrange( term.getNumResources()) if not any( term.getResourceURI( idx).startswith(s) for s in skip_uris)]
                if term.getBiologicalQualifierType() in [libsbml.BQB_IS, libsbml.BQB_IS_HOMOLOG_TO]:
                    bqbiol_is.extend( uris)
                elif term.getBiologicalQualifierType() == libsbml.BQB_HAS_PART:
                    bqbiol_has_part.extend( uris) 
                elif term.getBiologicalQualifierType() == libsbml.BQB_HAS_VERSION:
                    bqbiol_has_version.extend( uris) 
        
        sbo = species.getSBOTerm()
        if sbo == -1:
            sbo = None;
            sbo_str = None;
        else:
            sbo_str = "SBO:{0:07d}".format( sbo)
        
        annotation = {}
        for prefix in PREFIXES:
            annotation[ prefix.strip()] = species.getName().count( prefix)

        if species.getCompartment() == "default":
            compartment = None
            compartment_id = None
        else:
            compartment = model.getCompartment( species.getCompartment()).getName().lower().strip()
            compartment_id = species.getCompartment() 
        
        node_data = { "type" : "species",
                       "id" : prefix + species.getId(),
                       "name" : species.getName(),
                       "compartment" : compartment,
                       "compartment_id" : compartment_id,
                       "bqbiol_is" : tuple( sorted( set( bqbiol_is))),
                       "bqbiol_has_part" : tuple( sorted( set( bqbiol_has_part))),
                       "bqbiol_has_version" : tuple( sorted( set( bqbiol_has_version))),
                       "sbo" : sbo,
                       "sbo_str" : sbo_str,
                       "participants" : [],
                       "participant_ids" : [],
                       "annotation" : annotation};
        
        graph.add_node( prefix + species.getId(), node_data)
                       
    # add species reactions
    for reaction in model.getListOfReactions():
        bqbiol_is = []
        bqbiol_has_part = []
        bqbiol_has_version = []
        if reaction.getCVTerms() != None:
            for term in reaction.getCVTerms():
                uris = [ term.getResourceURI( idx) for idx in xrange( term.getNumResources()) if not any( term.getResourceURI( idx).startswith(s) for s in skip_uris)]
                if term.getBiologicalQualifierType() in [libsbml.BQB_IS, libsbml.BQB_IS_HOMOLOG_TO]:
                    bqbiol_is.extend( uris)
                elif term.getBiologicalQualifierType() == libsbml.BQB_HAS_PART:
                    bqbiol_has_part.extend( uris) 
                elif term.getBiologicalQualifierType() == libsbml.BQB_HAS_VERSION:
                    bqbiol_has_version.extend( uris) 
                    
        sbo = reaction.getSBOTerm()
        if sbo == -1:
            sbo = None;
            sbo_str = None;
        else:
            sbo_str = "SBO:{0:07d}".format( sbo)
            bqbiol_is.append( "urn:miriam:SBO:{0:07d}".format( sbo))
        graph.add_node( prefix + reaction.getId(), 
                       { "type" : "reaction",
                       "id" : prefix + reaction.getId(),
                       "local_id" : reaction.getId(),
                       "name" : reaction.getName(), 
                       "compartment" : reaction.getCompartment(),
                       "bqbiol_is" : tuple( sorted( set( bqbiol_is))),
                       "bqbiol_has_part" : tuple( sorted( set( bqbiol_has_part))),
                       "bqbiol_has_version" : tuple( sorted( set( bqbiol_has_version))),
                       "sbo" : sbo,
                       "sbo_str" : sbo_str} )
        # add edges
        for i in xrange( model.getNumReactions()):
            reaction = model.getReaction(i);
            for r in xrange( reaction.getNumReactants()):
                graph.add_edge( prefix + reaction.getId(), prefix + reaction.getReactant(r).getSpecies(), type = "reactant")
            for p in xrange( reaction.getNumProducts()):
                graph.add_edge( prefix + reaction.getId(), prefix + reaction.getProduct(p).getSpecies(), type = "product")
            for m in xrange( reaction.getNumModifiers()):
                graph.add_edge( prefix + reaction.getId(), prefix + reaction.getModifier(m).getSpecies(), type = "modifier")
    
    if ignore_participant_graph:
        return graph
    else:
        # participant graph
        participant_graph = networkx.DiGraph()
        graph_w_participant_edges = graph.copy()
        # add participant links
        for i in xrange( model.getNumSpecies()):
            species = model.getSpecies(i);
            graph_node = graph.node[ prefix + species.getId()]
            for participant in get_participants_species( species, prefix, model):
                # add participant graph edge
                participant_graph.add_edge( prefix + species.getId(), prefix + participant, type = "participant")
                graph_w_participant_edges.add_edge( prefix + species.getId(), prefix + participant, type = "participant")
                # add participant node information to 
                graph_node["participant_ids"].append( prefix + participant)
                graph_node["participants"].append( graph.node[prefix + participant])
                graph_node["bqbiol_has_part"] = tuple( sorted( set( list( graph.node[prefix + participant]["bqbiol_has_part"]) + list( graph_node["bqbiol_has_part"]))))

        return graph, participant_graph, graph_w_participant_edges

########################################################################
    
def bqbiol_is_map( graph):
    "returns a dictionary mapping of uri to node ids {uri : set( node ids)}"
    signature_map = {}
    for i in graph.nodes():
        node = graph.node[i]
        if signature_map.get( node["bqbiol_is"]) == None:
            signature_map[node["bqbiol_is"]] = [i]
        else:
            signature_map[node["bqbiol_is"]].append( i)
    return signature_map

def get_all_bqbiol_is_uris( graph):
    """ Returns all bqbiol_is uris from a graph """
    unique_ids = set()
    for n in graph.nodes( data = True):
        if n[1].get("bqbiol_is") and n[1].get("bqbiol_is") != ():
            unique_ids.update( n[1].get("bqbiol_is"))    
    return unique_ids

########################################################################
########################################################################

def find_nodes( graph, attribute, value):
    return [ n[1] for n in graph.nodes( data = True ) if n[1].get( attribute) != None and n[1][attribute] == value]
   
########################################################################
########################################################################

def filter_graph_remove_species_wo_bqbiol_is( graph):
    "Remove species without bqbiol_is"
    graph_cpy = graph.copy()
    remove_nodes = []
    for node in graph_cpy.nodes( data = True):
        # require nothing for reaction
        n = node[1]
        if n["type"] != "reaction" and n["bqbiol_is"] == ():
            remove_nodes.append( node[0])
    graph_cpy.remove_nodes_from( remove_nodes)
    graph_cpy.name = graph.name + "-REMOVED-SPECIES-WO-BQBIOL-IS"
    graph_cpy.file_name = None
    return graph_cpy

def filter_graph_remove_isolated_nodes( graph):
    "Remove nodes without connections"
    graph_cpy = graph.copy()
    graph_cpy.remove_nodes_from( networkx.isolates( graph))
    graph_cpy.name = graph.name + "-NO-ISOLATED-NODES"
    graph_cpy.file_name = None
    return graph_cpy

def filter_graph_remove_isolated_participants( graph):
    """Remove nodes without connections that are participants - 
     keep isolated nodes that are not participatns"""
    graph_cpy = graph.copy()
    isolates = set( networkx.isolates( graph))
    participants = set( [p for parts in [ [p ["id"] for p in n[1].get("participants")] for n in graph.nodes( data = True) if n[1].get("participants")] for p in parts ])
    graph_cpy.remove_nodes_from( isolates.intersection( participants))
    graph_cpy.name = graph.name + "-NO-ISOLATED-NODES"
    graph_cpy.file_name = None
    return graph_cpy
    
def filter_graph_remove_reactions_wo_sbo( graph):
    "Remove reactions without bqbiol_is"
    graph_cpy = graph.copy()
    remove_nodes = []
    for node in graph_cpy.nodes( data = True):
        # require nothing for reaction
        n = node[1]
        if n["type"] != "species" and n["sbo"] == None:
            remove_nodes.append( node[0])
    graph_cpy.remove_nodes_from( remove_nodes)
    graph_cpy.name = graph.name + "-REMOVED-REACTIONS-WO-SBO"
    graph_cpy.file_name = None
    return graph_cpy
    
def filter_reactions( graph):
    "remove all nodes that are NOT a reaction"
    graph_cpy = graph.copy()
    non_reaction_ids = [ n[0] for n in graph_cpy.nodes( data = True) if n[1]["type"] != "reaction"]
    graph_cpy.remove_nodes_from( non_reaction_ids)
    graph_cpy.name = graph.name + "-REACTIONS"
    graph_cpy.file_name = None
    return graph_cpy

def filter_reactions_sbo( graph):
    "Remove all nodes that are NOT reactions with SBO"
    graph_cpy = graph.copy()
    non_reaction_ids = [ n[0] for n in graph_cpy.nodes( data = True) if n[1]["type"] != "reaction" or n[1]["sbo"] == None]
    graph_cpy.remove_nodes_from( non_reaction_ids)
    graph_cpy.name = graph.name + "-SBO-REACTIONS"
    graph_cpy.file_name = None
    return graph_cpy

def filter_species( graph):
    "Remove all nodes that are NOT species"
    graph_cpy = graph.copy()
    non_species_ids = [ n[0] for n in graph_cpy.nodes( data = True) if n[1]["type"] != "species"]
    graph_cpy.remove_nodes_from( non_species_ids)
    graph_cpy.name = graph.name + "-SPECIES"
    graph_cpy.file_name = None
    return graph_cpy

def filter_species_bqbiol_is( graph):
    "Remove all nodes that are NOT species with bqbiol_is"
    graph_cpy = graph.copy()
    non_bqbiol_is_species_ids = [ n[0] for n in graph_cpy.nodes( data = True) if n[1]["type"] != "species" or n[1]["bqbiol_is"] == ()]
    graph_cpy.remove_nodes_from( non_bqbiol_is_species_ids)
    graph_cpy.name = graph.name + "-BQBIOL-IS-SPECIES"
    graph_cpy.file_name = None
    return graph_cpy
    
def filter_species_complex( graph):
    "Removes all nodes that are not complex - don't have participants"
    graph_cpy = graph.copy()
    non_complexes = [n[0] for n in graph.nodes( data = True) if not n[1].get("participants")]
    graph_cpy.remove_nodes_from( non_complexes)
    graph_cpy.name = graph.name + "-COMPLEXES"
    graph_cpy.file_name = None
    return graph_cpy

def filter_species_complex2( graph):
    "REmoves all nodes that are not complex - do not have sbo == 253"
    graph_cpy = graph.copy()
    non_complexes = [n[0] for n in graph.nodes( data = True) if not n[1].get("sbo") or n[1]["sbo"] != 253]
    graph_cpy.remove_nodes_from( non_complexes)
    graph_cpy.name = graph.name + "-COMPLEXES"
    graph_cpy.file_name = None
    return graph_cpy

########################################################################
########################################################################

def run_analysis( graph, export_file = None):
    """ Collects some simple statistics about the graph """
    import pandas
    print("%s:%s: run_analysis" % (now(), graph.name))
    species = filter_species( graph)
    reactions = filter_reactions( graph)
    edges = [n[2] for n in graph.edges( data = True)]
    isolated_nodes = networkx.isolates( graph)
    
    
    print("%s:%s: Computing statistics" % (now(), graph.name))
    d = {"name" : graph.name,
    
         "# nodes" : len( graph.nodes()),
         "# species" : len( species.nodes()),
         "# reactions" : len( reactions.nodes()),
    
         "# edges" : len( edges),
         "# edges reactant" : len( [ e for e in edges if e["type"] == "reactant"]),
         "# edges product" :  len( [ e for e in edges if e["type"] == "product"]), 
         "# edges modifier" : len( [ e for e in edges if e["type"] == "modifier"]),
    
         "# compartments" : len(set( [species.node[s]["compartment_id"] for s in species.nodes() if species.node[s]["compartment_id"]])),
         "# unique compartment names" : len(set( [species.node[s]["compartment"] for s in species.nodes() if species.node[s]["compartment"]])),
        
         "# isolated nodes" : len(isolated_nodes),
         "# isolated subgraphs" : len( list( networkx.connected_component_subgraphs( graph)))}
    data = pandas.Series(d)
    print("%s:%s: Results" % (now(), graph.name))
    print( data)
    if export_file:
        print("%s:%s: Exporting %s" % (now(), graph.name, export_file))
        data.to_pickle( export_file)
    
    print("%s:%s: Computing isolated nodes" % (now(), graph.name))
    isolates = set( networkx.isolates( graph))
    participants = set( [p for parts in [ [p ["id"] for p in n[1].get("participants")] for n in graph.nodes( data = True) if n[1].get("participants")] for p in parts ])
    real_isolates = isolates.difference( participants) # we have to discount those that are participants in a complex
    d["isolates # nodes"] = len( real_isolates)
    d["isolates # species"] = len( [n for n in real_isolates if graph.node[n]["type"] == "species"])
    d["isolates # reactions"] = len( [n for n in real_isolates if graph.node[n]["type"] == "reaction"])
    
    
    print("%s:%s: Computing subgraphs" % (now(), graph.name))
    # compute new graph with participant links
    participant_edges = []
    subgraphs = None
    for n1 in graph.nodes(data=True):
        if "participants" in n1[1] and n1[1]["participants"] != []:
            participant_edges.extend( [(n1[1]["id"], n2["id"]) for n2 in n1[1]["participants"]])
    if participant_edges != []:
        graph = graph.copy()
        [graph.add_edge( e[0], e[1], type = "participant") for e in participant_edges]
        subgraphs = list( networkx.connected_component_subgraphs( graph))
    elif subgraphs == None:
        subgraphs = list( networkx.connected_component_subgraphs( graph))
    nr_nodes = [ len( s.nodes()) for s in subgraphs]
    nr_edges = [ len( s.edges()) for s in subgraphs]
        
    d["subgraphs # subgraphs"] = len( subgraphs)
    
    d["subgraphs # nodes min"] = min(nr_nodes)
    d["subgraphs # nodes mean"] = numpy.mean( nr_nodes)
    d["subgraphs # nodes median"] = numpy.median( nr_nodes)
    d["subgraphs # nodes max"] = max( nr_nodes)
    d["subgraphs nodes histogram"] = collections.Counter( nr_nodes)
    
    d["subgraphs # edges min"] = min(nr_nodes)
    d["subgraphs # edges mean"] = numpy.mean( nr_nodes)
    d["subgraphs # edges median"] = numpy.median( nr_nodes)
    d["subgraphs # edges max"] = max( nr_nodes)
    d["subgraphs edges histogram"] = collections.Counter( nr_edges)
    
    data = pandas.Series(d)
    print("%s:%s: Results" % (now(), graph.name))
    print( data)
    if export_file:
        print("%s:%s: Exporting %s" % (now(), graph.name, export_file))
        data.to_pickle( export_file)
    return data

########################################################################
########################################################################

def run_analysis_signatures( graph, export_file = None, d = {}):
    """ Collects some statistics about the graphs names, bqbiol is signatures etc
    This takes a long time at this point. Use carefully"""
    print("%s:%s: run_analysis_signatures" % (now(), graph.name))
    import pandas
    
    species = filter_species( graph)
    reactions = filter_reactions( graph)
    
    
    if not "name" in d.keys():
        d["name"] = graph.name
    
    ## names
    print("%s:%s: Computing name statistics" % (now(), graph.name))
    species_names = [ species.node[n]["name"].lower() for n in species if species.node[n]["name"] != ""]
    d["species % have name"] = 100. * len(species_names) / len( species.nodes())
    d["species # unique names"] = len(set(species_names))
    species_clean_names = set([ clean_name(species.node[n]["name"]) for n in species if species.node[n]["name"] != ""])
    d["species # unique clean names"] = len(species_clean_names)
    species_clean_names2 = set([ clean_name2(species.node[n]["name"]) for n in species if species.node[n]["name"] != ""])
    d["species # unique clean names2"] = len(species_clean_names2)
    similar_names = []
    for name in species_clean_names2:
        similar = filter( lambda n: fuzzywuzzy.fuzz.ratio( name, n) > 90, species_clean_names2)
        similar_names.append( set( [name] + similar))
        similar_names = merge( similar_names)
    d["species # similar unique clean names2"] = len( similar_names)
    
    print("%s:%s: Computing bqbiol_is statistics species" % (now(), graph.name))
    species_bqbiolis = [ species.node[n]["bqbiol_is"] for n in species if species.node[n]["bqbiol_is"]]
    species_bqbiolis_signature_unique = set( species_bqbiolis)
    species_bqbiolis_terms = set( [ b for n in species for b in species.node[n]["bqbiol_is"]])
    
    d["species % have bqbiol_is"] = 100. * len( species_bqbiolis) / float( len(species))
    d["species # unique bqbiol_is signatures"] = len( species_bqbiolis_signature_unique)
    d["species # unique bqbiol_is terms"] = len( species_bqbiolis_terms)
    
    species_bqbiol_has_part = [ species.node[n]["bqbiol_has_part"] for n in species if species.node[n]["bqbiol_has_part"]]
    species_bqbiol_has_part_signature_unique = set( species_bqbiol_has_part)
    species_bqbiol_has_part_terms = set( [ b for n in species for b in species.node[n]["bqbiol_has_part"]])
    
    d["species % have bqbiol_has_part"] = 100* len( species_bqbiol_has_part) / float( len(species))
    d["species # unique bqbiol_has_part signatures"] = len( species_bqbiol_has_part_signature_unique)
    d["species # unique bqbiol_has_part terms"] = len( species_bqbiol_has_part_terms)
    
    print("%s:%s: Computing bqbiol_is statistics reactions" % (now(), graph.name))
    reactions_uri = [ reactions.node[n]["bqbiol_is"] for n in reactions if reactions.node[n]["bqbiol_is"]]
    reactions_uri_signature_unique = set( reactions_uri)
    reactions_bqbiol_terms = [ b for n in reactions for b in reactions.node[n]["bqbiol_is"]]
    reactions_bqbiol_terms_known = [ t for t in reactions_bqbiol_terms if sbo_go_name_known(t)]
    reactions_bqbiol_terms_set = set( reactions_bqbiol_terms)
    reactions_bqbiol_terms_known_set = set(reactions_bqbiol_terms_known)
    unknown_terms = reactions_bqbiol_terms_set.difference( reactions_bqbiol_terms_known_set)
    d["reactions % have bqbiol_is"] = 100* len( reactions_uri) / float( len(reactions))
    d["reactions # unique bqbiol_is signatures"] = len( reactions_uri_signature_unique)
    d["reactions # unique bqbiol_is terms"] = len(reactions_bqbiol_terms_set)
    d["reactions # unique bqbiol_is terms known SBO/GO terms"] = len(reactions_bqbiol_terms_set)
    d["reactions # unique bqbiol_is terms unknown SBO/GO terms"] = len( unknown_terms)
    d["reactions bqbiol_is terms histogram"] = collections.Counter( reactions_bqbiol_terms_known)

    data = pandas.Series(d)
    print("%s:%s: Results" % (now(), graph.name))
    print( data)
    if export_file:
        print("%s:%s: Exporting %s" % (now(), graph.name, export_file))
        data.to_pickle( export_file)

    
########################################################################
########################################################################

def run_analysis_isolated_nodes( graph):
    print( "\n\nrun_analysis_isolated_nodes(%s)" % graph.name)
    isolates = set( networkx.isolates( graph))
    participants = set( [p for parts in [ [p ["id"] for p in n[1].get("participants")] for n in graph.nodes( data = True) if n[1].get("participants")] for p in parts ])
    
    real_isolates = isolates.difference( participants) # we have to discount those that are participants in a complex

    print( "{} isolated nodes (ignoring participant nodes) ({} isolated species, {} isolated reactions)".format(
        len( real_isolates),
        len( [n for n in real_isolates if graph.node[n]["type"] == "species"]),
        len( [n for n in real_isolates if graph.node[n]["type"] == "reaction"])))
    print( "{} isolated nodes (including isolated participant nodes) ({} isolated species, {} isolated reactions)".format(
        len( isolates),
        len( [n for n in isolates if graph.node[n]["type"] == "species"]),
        len( [n for n in isolates if graph.node[n]["type"] == "reaction"])))

########################################################################
########################################################################

def run_analysis_subgraphs( graph, subgraphs = None):
    """ Compute some statistics for subgraphs: min,max,median """
    print( "run_analysis_subgraphs( %s)" % graph.name)
    # compute new graph with participant links
    participant_edges = []
    for n1 in graph.nodes(data=True):
        if "participants" in n1[1] and n1[1]["participants"] != []:
            participant_edges.extend( [(n1[1]["id"], n2["id"]) for n2 in n1[1]["participants"]])
    if participant_edges != []:
        graph = graph.copy()
        [graph.add_edge( e[0], e[1], type = "participant") for e in participant_edges]
        subgraphs = list( networkx.connected_component_subgraphs( graph))
    elif subgraphs == None:
        subgraphs = list( networkx.connected_component_subgraphs( graph))
        

    nr_nodes = [ len( s.nodes()) for s in subgraphs]
    nr_edges = [ len( s.edges()) for s in subgraphs]
    
    print( "{} # subgraphs".format( len( subgraphs)))
    print( "{}/{}/{}/{} min/mean/median/max # nodes per subgraph".format( min(nr_nodes), numpy.mean( nr_nodes), numpy.median( nr_nodes), max( nr_nodes)))
    print( "{}/{}/{}/{} min/mean/median/max # edges per subgraph".format( min(nr_edges), numpy.mean( nr_edges), numpy.median( nr_nodes), max( nr_edges)))
    print()
    print( "# nodes per subgraph statistics: {}".format( collections.Counter( nr_nodes)))
    print( "# edges per subgraph statistics: {}".format( collections.Counter( nr_edges)))


    subgraphs_no_isolates = [ s for s in subgraphs if len(s.nodes()) > 1]
    nr_nodes_subgraphs_no_isolates = [ len( s.nodes()) for s in subgraphs_no_isolates]
    nr_edges_subgraphs_no_isolates = [ len( s.edges()) for s in subgraphs_no_isolates]
    print( "\n--\n")
    print( "{} # subgraphs no isolated nodes".format( len( subgraphs_no_isolates)))
    print( "{}/{}/{}/{} min/mean/median/max # nodes per subgraphs no isolated nodes".format( min( nr_nodes_subgraphs_no_isolates), numpy.mean( nr_nodes_subgraphs_no_isolates), numpy.median( nr_nodes_subgraphs_no_isolates), max( nr_nodes_subgraphs_no_isolates)))
    print( "{}/{}/{}/{} min/mean/median/max # edges per subgraphs no isolated nodes".format( min( nr_edges_subgraphs_no_isolates), numpy.mean( nr_edges_subgraphs_no_isolates), numpy.median( nr_edges_subgraphs_no_isolates), max( nr_edges_subgraphs_no_isolates)))
    print()
    print( "# nodes per subgraph (no isolated nodes) statistics: {}".format( collections.Counter( nr_nodes_subgraphs_no_isolates)))
    print( "# edges per subgraph (no isolated nodes) statistics: {}".format( collections.Counter( nr_edges_subgraphs_no_isolates)))
        
########################################################################
########################################################################
        
def run_analysis_complex_participants( graph, participant_graph):
    
    node_dict = { n[0] : n[1] for n in graph.nodes(data = True) }
    edges_dict = { n: [] for n in node_dict.keys()}
    for e in graph.edges( data = True):
        edges_dict[e[0]].append( (e[1],e[2]["type"]))
        edges_dict[e[1]].append( (e[0],e[2]["type"]))
    
    reaction_participants = set( [ p for e in graph.edges() for p in e])

    # sbo complexes
    complexes = [n[1] for n in graph.nodes() if n[1]["type"] == "species" and n[1]["sbo"] == 253]     
    complexes_ids = set( [ c["id"] for c in complexes])
    assert( len( complexes) == len( complexes_ids))
    print( "{} total # of complexes (sbo == 253)".format( len( complexes)))    
    
    # complexes based on Participant edge
    complexes2 = set( [ e[0] for e in participant_graph.edges()])
    complexes2_participant = set( [ e[1] for e in participant_graph.edges()]) # participants of complexes
    print( "{} total # of complexes (complex in a complex relationship with some participant)".format( len( complexes2)))    
    print( "{} total # of unique participants".format( len( complexes2_participant)))
    
    # complexes part of reaction
    complexes_in_reaction = complexes_ids.intersection( reaction_participants)
    complexes_not_in_reaction = complexes_ids.difference( reaction_participants)
    print( "{}/{} of complexes are part of a reaction ({}/{} are not)".format( 
        len( complexes_in_reaction), 
        len( complexes_ids),    
        len( complexes_not_in_reaction), 
        len( complexes_ids)))    
    
    # participants part of reaction
    complexes_participant_in_reaction = complexes2_participant.intersection( reaction_participants)
    complexes_participant_not_in_reaction = complexes2_participant.difference( reaction_participants)
    print( "{}/{} of participants are part of a reaction ({}/{} are not)".format( 
        len( complexes_participant_in_reaction), 
        len( complexes2_participant),    
        len( complexes_participant_not_in_reaction), 
        len( complexes2_participant)))    
        
    complexes_participants_in_other_complexes = complexes_ids.intersection( complexes2_participant)
    print( "{} complexes participate in other complexes".format( len( complexes_participants_in_other_complexes))) 
    
    multiple_complex_edge_participant = [n for n, c in collections.Counter( [ e[1] for e in participant_graph.edges()]).items() if c > 1]
    print( "{} participants participate in multiple complexes".format( len(multiple_complex_edge_participant))) 
    
    ## some annotation information
    complexes_wo_bqbiol_is = [  c for c in complexes_ids if graph.node[c]["bqbiol_is"] == ()]
    print( "{}/{} complexes w/o bqbiol_is".format( len( complexes_wo_bqbiol_is), len( complexes_ids))) 
    participants_wo_bqbiol_is = [  p for p in complexes2_participant if graph.node[p]["bqbiol_is"] == ()]
    print( "{}/{} participants w/o bqbiol_is".format( len( participants_wo_bqbiol_is), len( complexes2_participant))) 

########################################################################
########################################################################

def precision_recall_f_score( tp, fp, fn):
    if len( tp) == 0 and len( fp) == 0:
        precision = 0
    else:
        precision = len( tp) / float( len( tp) + len( fp))
    if len( tp) == 0 and len( fn) == 0:
        recall = 0
    else:
        recall = len( tp) / float( len( tp) + len( fn))
    if precision == 0 and recall == 0:
        f_score =  0.0
    else:
        f_score =  2.0 * (precision * recall) / (precision + recall)
    return precision, recall, f_score

########################################################################
########################################################################

def set_overlap( set_1, set_2, equal_fn):
    r_1 = set()
    r_2 = set()
    for e1 in set_1:
        e2s = filter( lambda e2: equal_fn( e1, e2), set_2)
        if e2s:
            r_2 = r_2.union( e2s)
            r_1.add( e1)
    return r_1, r_2

def list_overlap( list_1, list_2, equal_fn):
    """ Returns indices of overlapping elements"""
    indices_1 = set()
    indices_2 = set()
    for i_1, e1 in enumerate( list_1):
        is_2 = [i for i, e2 in enumerate(list_2) if equal_fn( e1, e2)]
        if is_2 != []:
            indices_2.update( is_2)
            indices_1.add( i_1)
    return indices_1, indices_2

def list_intersect( list_1, list_2):
    l_1 = list_1[:]
    l_2 = list_2[:]
    result = []
    
    while len(l_1) > 0:
        e1 = l_1.pop()
        try:
            idx = l_2.index( e1)
        except:
            idx = None
        if idx != None:
            l_2.remove( e1)
            result.append( e1)
    return result
    
assert( list_intersect([1,2,3],[4,5]) == [])
assert( list_intersect([1,2,3],[1,5]) == [1])
assert( list_intersect([1,2,3,1],[1,5]) == [1])
assert( list_intersect([1,2,3,1],[1,1]) == [1,1])

def list_difference( list_1, list_2):
    l_1 = list_1[:]
    for e2 in list_2:
        try:
            l_1.remove(e2)
        except:
            pass
        
    return l_1
    
assert( list_difference([1,2,3,1],[5,6]) == [1,2,3,1])
assert( list_difference([1,2,3,1],[1,6]) == [2,3,1])
assert( list_difference([1,2,3,1],[1,1,6]) == [2,3])

def list_find( el, listt, equal_fn):
    for el2 in listt:
        if equal_fn( el, el2):
            return el2
    return None
    
def list_difference2( list_1, list_2, equal_fn):
    "returns those elements of list_1 which are not in list_2 according to equal_fn"
    result = []
    for e in list_1:
        if not list_find( e, list_2, equal_fn):
            result.append( e)
    return result


def list_reduce2( list_1, equal_fn):
    result = []
    elements_remaining = list_1[:]
    while elements_remaining:
        el = elements_remaining.pop()
        result.append( el)
        new_elements_remaining = []
        for el2 in elements_remaining:
            if not equal_fn( el, el2):
                new_elements_remaining.append( el2)
        elements_remaining = new_elements_remaining
    return result

assert( list_reduce2([1,"1",2,"2"], lambda e1, e2: str( e1) == str( e2)) == ['2','1'])

def merge( sets):
  "merges sets which are disjoint"
  merged = 1
  while merged:
    merged = 0
    results = []
    while sets:
      common, rest = sets[0], sets[1:]
      sets = []
      for x in rest:
        if x.isdisjoint(common):
          sets.append(x)
        else:
          merged = 1
          common |= x
      results.append(common)
    sets = results
  return sets
    
def analyse_set_overlap( set_1, set_2, equal_fn = operator.eq):
    res_1, res_2 = set_overlap( set_1, set_2, equal_fn)
    if len( set_2) == 0:
        precision = 0
    else:
        precision = 100.0 * len( res_2) / float( len( set_2))
    if len( set_1) == 0:
        recall = 0
    else:
        recall = 100.0 * len( res_1) / float( len( set_1))
    if precision == 0 and recall == 0:
        f_score =  0.0
    else:
        f_score =  2.0 * (precision * recall) / (precision + recall)
    return res_1, res_2, precision, recall, f_score
    
def analyse_list_overlap( list_1, list_2, equal_fn = operator.eq):
    res_1, res_2 = list_overlap( list_1, list_2, equal_fn)
    if len( list_2) == 0:
        precision = 0
    else:
        precision = 100.0 * len( res_2) / float( len( list_2))
    if len( list_1) == 0:
        recall = 0
    else:
        recall = 100.0 * len( res_1) / float( len( list_1))
    if precision == 0 and recall == 0:
        f_score =  0.0
    else:
        f_score =  2.0 * (precision * recall) / (precision + recall)
    return res_1, res_2, precision, recall, f_score

def tuple_eq_empty_not_eq( t_1, t_2):
    """ those which are empty are in fact not equal"""
    return len( t_1) > 0 and  t_1 == t_2
        
def tuple_overlaps( t_1, t_2):
    return len( set(t_1).intersection( t_2)) > 0

def tuple_overlaps_sbo_is_a( t_1, t_2):
    if tuple_overlaps( t_1, t_2):
        return True
    else:
        sbo_terms_1 = get_sbo_terms( t_1)
        sbo_terms_2 = get_sbo_terms( t_2)
        for s1 in sbo_terms_1:
            for s2 in sbo_terms_2:
                if sbo_is_a2( s1, s2):
                    return True
                    
def name_approx_equal( n1, n2):
    return fuzzywuzzy.fuzz.ratio( n1, n2) > 90

########################################################################
########################################################################

def nm_name_equal( n1, n2):
    "Checks if name is the same"
    return n1["name"].lower() == n2["name"].lower()

def nm_name_equal_w_participants( n1, n2):
    "Checks if name and names of participants overlap"
    names_1 = [n1["name"].lower()] + [ p["name"].lower() for p in n1["participants"]]
    names_2 = [n2["name"].lower()] + [ p["name"].lower() for p in n2["participants"]]
    return len( set( names_1).intersection( names_2)) > 0

def nm_name_clean_equal( n1, n2):
    "Checks if clean name is the same"
    return remove_prefixes( n1["name"].lower()) == remove_prefixes( n2["name"].lower())

def nm_name_clean_equal_w_participants( n1, n2):
    "Checks if name and names of participants overlap"
    clean_names_1 = [remove_prefixes( n1["name"].lower())] + [ remove_prefixes( p["name"].lower()) for p in n1["participants"]]
    clean_names_2 = [remove_prefixes( n2["name"].lower())] + [ remove_prefixes( p["name"].lower()) for p in n2["participants"]]
    return len( set( clean_names_1).intersection( clean_names_2)) > 0

def nm_name_clean2_equal( n1, n2):
    "Checks if clean name is the same"
    return clean_name2( n1["name"]) == clean_name2( n2["name"])
    
def nm_name_clean_approx( n1, n2):
    return fuzzywuzzy.fuzz.ratio( clean_name2( n1["name"]), clean_name2( n2["name"])) > 90
    
def nm_name_clean_approx_w_participants( n1, n2):
    clean_names_1 = [ re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( n1["name"].lower()))] + [ re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( p["name"].lower())) for p in n1["participants"]]
    clean_names_2 = [ re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( n2["name"].lower()))] + [ re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( p["name"].lower())) for p in n2["participants"]]
    for name_1 in clean_names_1:
        if list_find( name_1, clean_names_2, lambda name_1, name_2: fuzzywuzzy.fuzz.ratio( name_1, name_2) > 90):
            return True
    return False
    
def nm_gene_id_intersect( n1, n2):
    set_1 = map_gene_simstring( n1["name"])
    set_2 = map_gene_simstring( n2["name"])
    return set_1 and set_2 and len( set_1.intersection( set_2)) > 0

def nm_gene_id_intersect_w_participants( n1, n2):
    sets_1 = filter( bool, [map_gene_simstring(n) for n in [ n1["name"]] + [ p["name"] for p in n1["participants"]]])
    sets_2 = filter( bool, [map_gene_simstring(n) for n in [ n2["name"]] + [ p["name"] for p in n2["participants"]]])
    for s1 in sets_1:
        for s2 in sets_2:
            if len( s1.intersection( s2)) > 0:
                return True
    return False
    
def nm_name_clean_approx_OR_gene_id_intersect( n1, n2):
    return nm_name_clean_approx( n1, n2) or nm_gene_id_intersect( n1, n2)

def nm_name_clean_approx_OR_gene_id_intersect_w_participants( n1, n2):
    return nm_name_clean_approx_w_participants( n1, n2) or nm_gene_id_intersect_w_participants( n1, n2)

def nm_bqbiol_is_equal( n1, n2):
    "Checks if the bqbiol_is are the same"
    return n1["bqbiol_is"] and n2["bqbiol_is"] and n1["bqbiol_is"] == n2["bqbiol_is"]

def nm_bqbiol_is_equal_w_participants( n1, n2):
    "Checks if the bqbiol_is are the same - also checks participants"
    sets_1 = filter( bool, [set(n1["bqbiol_is"])] +  [ set(p["bqbiol_is"]) for p in n1["participants"]])
    sets_2 = filter( bool, [n2["bqbiol_is"]] +  [ p["bqbiol_is"] for p in n2["participants"]])
    for s1 in sets_1:
        for s2 in sets_2:
            if len( s1.intersection( s2)) > 0:
                return True
    return False

def nm_bqbiol_is_overlaps( n1, n2):
    "Checks if the bqbiol_is are the same"
    return n1["bqbiol_is"] and n2["bqbiol_is"] and len( set( n1["bqbiol_is"]).intersection(  set( n2["bqbiol_is"]))) > 0

def nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    "Checks if the bqbiol_is are the same"
    if nm_bqbiol_is_overlaps( n1, n2):
        return True
    elif n1["bqbiol_is"] and n2["bqbiol_is"]:
        sbo_terms_1 = get_sbo_terms( n1["bqbiol_is"])
        sbo_terms_2 = get_sbo_terms( n2["bqbiol_is"])
        for s1 in sbo_terms_1:
            for s2 in sbo_terms_2:
                if sbo_is_a2( s1, s2):
                    return True
    return False

def nm_bqbiol_is_overlaps_w_participants( n1, n2):
    "Checks if the bqbiol_is overlaps - also checks participants"
    set_1 = set( n1["bqbiol_is"])
    if n1.get("participants"):
        [set_1.update( p["bqbiol_is"]) for p in n1["participants"]]
    set_2 = set( n2["bqbiol_is"])
    if n2.get("participants"):
        [set_2.update( p["bqbiol_is"]) for p in n2["participants"]]
    if len( set_1.intersection( set_2)) > 0:
        return True
    else:
        return False

def nm_bqbiol_is_has_part_overlaps( n1, n2):
    "Checks if the bqbiol_is and bqbiol_has_part overlaps"
    uris_1 = set()
    if n1["bqbiol_is"]:
        uris_1.update( n1["bqbiol_is"])
    if n1["bqbiol_has_part"]:
        uris_1.update( n1["bqbiol_has_part"])
    uris_2 = set()
    if n2["bqbiol_is"]:
        uris_2.update( n2["bqbiol_is"])
    if n1["bqbiol_has_part"]:
        uris_2.update( n2["bqbiol_has_part"])
    return len( uris_1.intersection(  uris_2)) > 0

def nm_sbo_equal( n1, n2):
    "Only works on reactions"
    return n1["sbo"] and n2["sbo"] and n1["sbo"] == n2["sbo"]

def nm_sbo_is_a( n1, n2):
    "Only works on reactions"
    return n1["sbo_str"] and n2["sbo_str"] and sbo_is_a2( n1["sbo_str"], n2["sbo_str"])

################### name_clean + various reactions matches

def nm_name_clean_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_equal( n1, n2):
            return True
    else:
        return False

def nm_name_clean_w_participants_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
        return True
    elif n1["type"] == "species" and nm_name_clean_equal_w_participants( n1, n2):
        return True
    else:
        return False
        
def nm_name_clean_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_equal( n1, n2):
            return True
    else:
        return False

def nm_name_clean_w_participants_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_equal_w_participants( n1, n2):
            return True
    else:
        return False

def nm_name_clean_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_equal( n1, n2):
            return True
    else:
        return False

def nm_name_clean_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_equal_w_participants( n1, n2):
            return True
    else:
        return False

################### name_clean_approx + various reactions matches

def nm_name_clean_approx_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_approx( n1, n2):
            return True
    else:
        return False

def nm_name_clean_approx_w_participants_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
        return True
    elif n1["type"] == "species" and nm_name_clean_approx_w_participants( n1, n2):
        return True
    else:
        return False
        
def nm_name_clean_approx_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_approx( n1, n2):
            return True
    else:
        return False

def nm_name_clean_approx_w_participants_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_approx_w_participants( n1, n2):
            return True
    else:
        return False

def nm_name_clean_approx_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_approx( n1, n2):
            return True
    else:
        return False

def nm_name_clean_approx_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species" and nm_name_clean_approx_w_participants( n1, n2):
            return True
    else:
        return False

################### name_clean_approx or bqbiol_is_equal various reactions matches

def nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_equal( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_equal_w_participants_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx_w_participants( n1, n2) or nm_bqbiol_is_equal_w_participants( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_equal( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_equal_w_participants_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx_w_participants( n1, n2) or nm_bqbiol_is_equal_w_participants( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_equal( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_equal_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx_w_participants( n1, n2) or nm_bqbiol_is_equal_w_participants( n1, n2)):
            return True
    else:
        return False

################### name_clean_approx or bqbiol_is_overlaps various reactions matches

def nm_name_clean_approx_OR_bqbiol_is_overlaps_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_overlaps( n1, n2)):
            return True
    else:
        return False
        
def nm_name_clean_approx_OR_bqbiol_is_overlaps_w_participants_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx_w_participants( n1, n2) or nm_bqbiol_is_overlaps_w_participants( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_overlaps_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_overlaps( n1, n2)):
            return True
    else:
        return False
        
def nm_name_clean_approx_OR_bqbiol_is_overlaps_w_participants_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx_w_participants( n1, n2) or nm_bqbiol_is_overlaps_w_participants( n1, n2)):
            return True
    else:
        return False
        
def nm_name_clean_approx_OR_bqbiol_is_overlaps_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_overlaps( n1, n2)):
            return True
    else:
        return False
        
def nm_name_clean_approx_OR_bqbiol_is_overlaps_w_participants_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx_w_participants( n1, n2) or nm_bqbiol_is_overlaps_w_participants( n1, n2)):
            return True
    else:
        return False

################### name_clean_approx or bqbiol_is_overlaps various reactions matches

def nm_name_clean_approx_OR_bqbiol_is_bqbiol_is_has_parts_overlaps_AND_nm_bqbiol_is_equal( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_equal( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_has_part_overlaps( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_bqbiol_is_has_parts_overlaps_AND_nm_bqbiol_is_overlaps( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_has_part_overlaps( n1, n2)):
            return True
    else:
        return False

def nm_name_clean_approx_OR_bqbiol_is_bqbiol_is_has_parts_overlaps_AND_nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
    if n1["type"] != n2["type"]:
        return False
    elif n1["type"] == "reaction" and nm_bqbiol_is_overlaps_sbo_is_a( n1, n2):
            return True
    elif n1["type"] == "species"and (nm_name_clean_approx( n1, n2) or nm_bqbiol_is_has_part_overlaps( n1, n2)):
            return True
    else:
        return False

################### edge match exact
        
def edge_match_exact( e1, e2):
    "only edges"
    return e1["type"] ==  e2["type"]  

########################################################################
########################################################################
# nodes overlap max

def compute_nodes_overlap_max( graph_1, graph_2, node_match):
    """ computes a nodes in graph_2 matching with nodes in graph_1 according
        to node_match
        Returns - a dictionary of nodes  """
    nodes_2 = [ filter( lambda n2: node_match( graph_1.node[n1], graph_2.node[n2]), graph_2.nodes()) for n1 in graph_1.nodes()]
    return { n1: n2 for n1, n2 in zip( graph_1.nodes(), nodes_2) if n2 }

def get_nodes_overlap_max_result_precision_recall_f_score( graph_1, graph_2, matches):
    if len( graph_2) == 0:
        precision = 0;
    else:
        precision = len( set( itertools.chain(*matches.values()))) / float( len( graph_2))
    if len( graph_1) == 0:
        recall = 0
    else:
        recall = len( matches.keys()) / float( len( graph_1))
    if precision == 0 and recall == 0:
        f_score =  0.0
    else:
        f_score =  2.0 * (precision * recall) / (precision + recall)

    return 100.0 * precision, 100.0 * recall, 100.0 * f_score

def print_node_match_result( graph_1, graph_2, matches, node_match_name = "", export_matches = None):
    # print results
    precision, recall, f_score = get_nodes_overlap_max_result_precision_recall_f_score( graph_1, graph_2, matches)
    print( "{}: {:.2f} & {:.2f} & {:.2f} node overlap (precision/recall/f-score)".format( 
            node_match_name, precision, recall, f_score))
    # export text matches files
    if export_matches:
        with open( export_matches, "wt") as f:
            clean_names_map = { clean_name2( graph_1.node[k]["name"]) : k for k in matches.keys()}
            for n in sorted( clean_names_map.keys()):
                k = clean_names_map[n]
                if matches[k]:
                    f.write( "\n-------------------------------------------------------------\n")
                    f.write( n)
                    f.write( "\n--\n" )
                    names = set( [clean_name2( graph_2.node[v]["name"]) for v in matches[k]])
                    f.write( "\n".join(names))
        
def run_analysis_nodes_overlap_max( graph_1, graph_2, node_match, 
                                   export_results = False, 
                                   export_results_prefix = "results-nodes-overlap-max",
                                   ignore_existing = False):
    """ computes nodes overlap and prints statistics"""
    export_file = "%s__%s__%s__%s.pickle" % (export_results_prefix, graph_1.name, graph_2.name, node_match.__name__)
    if ignore_existing and os.path.exists( export_file):        
        print("%s:%s/%s:run_analysis_nodes_overlap_max:%s exists. using that one." % (now(),graph_1.name, graph_2.name, export_file))
        data = pickle.load( open( export_file, "rb"))
        graph_1, graph_2, matches = data[0], data[1], data[2]
    else:
        matches = compute_nodes_overlap_max( graph_1, graph_2, node_match)
    print_node_match_result( graph_1, graph_2, matches, node_match_name = node_match.__name__)
    if export_results and not( ignore_existing and os.path.exists( export_file)):
        print("%s:%s/%s:run_analysis_nodes_overlap_max:Exporting %s" % (now(),graph_1.name, graph_2.name, export_file))
        pickle.dump( [graph_1, graph_2, matches], open( export_file, "wb"))

def run_analyses_nodes_overlap_max( graph_1, 
                               graph_2, 
                               node_match_fns,
                               prefix = None,
                               n_jobs = None,
                               export_results = False,
                               export_results_prefix = "results-nodes-overlap-max"):
    """ computes nodes overlaps according to multiple node_match_fns and prints statistics """
    print( "-----------")
    print( "run_analyses_nodes_overlap_max %s/%s n_jobs=%s -- %s" % (graph_1.name, graph_2.name, n_jobs, node_match_fns))
        
    # compute the nodes of 2 that exist in 1 (ignoring edges)
    if n_jobs:
        with joblib.Parallel( n_jobs = n_jobs) as parallel:
            parallel( joblib.delayed( run_analysis_nodes_overlap_max) ( graph_1, graph_2, fn, export_results = export_results, export_results_prefix = export_results_prefix) 
            for fn in node_match_fns)        
    else:
        for nm in node_match_fns:
            run_analysis_nodes_overlap_max( graph_1, graph_2, nm, export_results = export_results, export_results_prefix = export_results_prefix)

########################################################################
########################################################################
# subgraph overlap max

def match_subgraph_max( graph, subgraph, node_match, edge_match = edge_match_exact, file_name = None):
    """ computes overlap for single subgraph"""                                                             
    assert( subgraph or file_name)
    if subgraph == None:
        subgraph = pickle.load( open( file_name, "rb"))
    graph_matcher = networkx.algorithms.isomorphism.GraphMatcher( graph, 
                                                                 subgraph, 
                                                                 node_match = node_match,
                                                                 edge_match = edge_match)
                                                            
    result = list( graph_matcher.subgraph_isomorphisms_iter())
    return result, subgraph

def match_graph_max( graph, file_name, node_match, edge_match = edge_match_exact):
    """ computes overlap for graph loaded from a file"""
    graph_2 = pickle.load( open( file_name, "rb"))
    subgraphs = list( networkx.connected_component_subgraphs( graph_2))
    graph_matchers = [networkx.algorithms.isomorphism.GraphMatcher( graph,
                                                                 subgraph, 
                                                                 node_match = node_match,
                                                                 edge_match = edge_match)
                                                                 for subgraph in subgraphs]
    
    results = [ (list( m.subgraph_isomorphisms_iter()), s) for m, s in zip( graph_matchers, subgraphs)]
    return results

def match_subgraphs_max( graph, 
                        subgraphs, 
                        node_match, 
                        edge_match = edge_match_exact, 
                        n_jobs = None,
                        file_names = None):
    """ computes overlap for subgraphs """                                                             
    # compute the nodes of 2 that exist in 1 (ignoring edges)
    
    assert( subgraphs or file_names)
    if file_names: # use the files instead of subgraphs if possible
        print( "Running match_subgraphs_max using file_names (individual graph files) n_jobs=%s" %(n_jobs))
        if n_jobs:
            with joblib.Parallel( n_jobs = n_jobs) as parallel:
                results = parallel( joblib.delayed( match_graph_max) ( graph, file_name, node_match = node_match, edge_match = edge_match) for file_name in file_names)
        else:
            results = [ match_graph_max( graph, file_name, node_match = node_match, edge_match = edge_match) for file_name in file_names]
        results = [r for result in results for r in result]
    else:
        print( "Running match_subgraphs_max using subgraphs n_jobs=%s" %(n_jobs))
        if n_jobs:
            with joblib.Parallel( n_jobs = n_jobs) as parallel:
                results = parallel( joblib.delayed( match_subgraph_max) ( graph, subgraph, node_match, edge_match) for subgraph in subgraphs)
        else:
            results = [ match_subgraph_max( graph, subgraph, node_match = node_match, edge_match = edge_match) for subgraph in subgraphs]
    results_matches = [r[0] for r in results]
    results_subgraphs = [r[1] for r in results]
    return results_matches, results_subgraphs

def subgraph_match_get_edges( subgraph, match, reverse_match, edge_signatures_1, edge_signatures_2):
    """ Computes matching edges from match_subgraphs results """
    m_edges = {}
    for e in subgraph.edges( data = True):
        # a bit of acrobatics to get around having to use digraph (which is buggy)
        signature_1_1 = (reverse_match[e[0]], reverse_match[e[1]], e[2]["type"])
        signature_1_2 = (reverse_match[e[1]], reverse_match[e[0]], e[2]["type"])
        signature_2_1 = (e[0], e[1], e[2]["type"])
        signature_2_2 = (e[1], e[0], e[2]["type"])
        assert signature_1_1 in edge_signatures_1 or signature_1_2 in edge_signatures_1
        assert not( signature_1_1 in edge_signatures_1 and signature_1_2 in edge_signatures_1)
        assert signature_2_1 in edge_signatures_2 or signature_2_2 in edge_signatures_2
        assert not( signature_2_1 in edge_signatures_2 and signature_2_2 in edge_signatures_2)
        if signature_1_1 in edge_signatures_1:
            signature_1 = signature_1_1
        else:
            signature_1 = signature_1_2
        if signature_2_1 in edge_signatures_2:
            signature_2 = signature_2_1
        else:
            signature_2 = signature_2_2
        
        m_edges[signature_1] = signature_2
        assert signature_1 in edge_signatures_1
        assert signature_2 in edge_signatures_2
    return m_edges

def compute_subgraphs_overlap_max( graph_1, graph_2, 
                                  node_match, 
                                  edge_match = edge_match_exact, 
                                  subgraphs_2 = None,
                                  n_jobs = None,
                                  export_results = False,
                                  export_results_prefix = "results-subgraphs-overlap-max",
                                  file_names = None,
                                  ignore_existing = False):
    """ compute the subgraphs in graph_1 isomorph to nodes in subgraphs of 2 """
    if export_results:
        export_file = "%s__%s__%s__%s__%s.pickle" % (export_results_prefix, graph_1.name, graph_2.name, node_match.__name__, edge_match.__name__)
    if export_results and ignore_existing and os.path.exists( export_file):
        print( "%s:%s/%s:compute_subgraphs_overlap_max:results exist %s, loading" % (now(), graph_1.name, graph_2.name, export_file))
        data = pickle.load( open( export_file, "rb"))
        graph_1, graph_2, results_subgraphs, results_matches = data[0], data[1], data[2], data[3]
        return results_matches, results_subgraphs
    
    if graph_2 and file_names == None and subgraphs_2 == None:
        subgraphs_2 = list( networkx.connected_component_subgraphs( graph_2))
    
    # Run!
    results_matches, results_subgraphs = match_subgraphs_max( graph_1, subgraphs_2, node_match = node_match, edge_match = edge_match, n_jobs = n_jobs, file_names = file_names)
    
    # export data
    if export_results:
        pickle.dump( [graph_1, graph_2, results_subgraphs, results_matches],
                    open( "%s__%s__%s__%s__%s.pickle" % (export_results_prefix, graph_1.name, graph_2.name, node_match.__name__, edge_match.__name__), "wb"))

    return results_matches, results_subgraphs

def get_subgraphs_overlap_max_results( graph_1, graph_2, results_subgraphs, results_matches, 
                                       species_1 = None, species_2 = None, reactions_1 = None, reactions_2 = None):
    """ takes results from matching and computes matches for nodes, edges, species, and reactions """
    if species_1 == None:
        species_1 = set(filter_species( graph_1).nodes())
    if species_2 == None:
        species_2 = set(filter_species( graph_2).nodes())
    if reactions_1 == None:
        reactions_1 = set( filter_reactions( graph_1).nodes())
    if reactions_2 == None:
        reactions_2 = set( filter_reactions( graph_2).nodes())
    

    # collet results for analysis
    matches_nodes_1 = set()
    matches_nodes_2 = set()    
    matches_edges_1 = set()
    matches_edges_2 = set()
    
    edge_signatures_1 = edge_signatures( graph_1)
    edge_signatures_2 = edge_signatures( graph_2)

    for subgraph_2, matches in zip( results_subgraphs, results_matches):
        for m in matches:
            matches_nodes_1 = matches_nodes_1.union( m.keys())
            matches_nodes_2 = matches_nodes_2.union( m.values())
            
            reverse_m = { v: k for k, v in m.iteritems()}
            m_edges = subgraph_match_get_edges( subgraph_2, m , reverse_m, edge_signatures_1, edge_signatures_2)            
            matches_edges_1 = matches_edges_1.union( m_edges.keys())
            matches_edges_2 = matches_edges_2.union( m_edges.values())

    species_1_matches = species_1.intersection( matches_nodes_1)
    species_2_matches = species_2.intersection( matches_nodes_2)
    reactions_1_matches = reactions_1.intersection( matches_nodes_1)
    reactions_2_matches = reactions_2.intersection( matches_nodes_2)

    return matches_nodes_1, matches_nodes_2, matches_edges_1, matches_edges_2, species_1_matches, species_2_matches, reactions_1_matches, reactions_2_matches

def get_subgraphs_overlap_max_results_precision_recall_f_score(graph_1, graph_2, results_subgraphs, results_matches, 
                                       species_1 = None, species_2 = None, reactions_1 = None, reactions_2 = None):
    """ Returns precision recall for nodes, species, reactions, edges as a dict """
    if species_1 == None:
        species_1 = set(filter_species( graph_1).nodes())
    if species_2 == None:
        species_2 = set(filter_species( graph_2).nodes())
    if reactions_1 == None:
        reactions_1 = set( filter_reactions( graph_1).nodes())
    if reactions_2 == None:
        reactions_2 = set( filter_reactions( graph_2).nodes())

    matches_nodes_1, matches_nodes_2, matches_edges_1, matches_edges_2, species_1_matches, species_2_matches, reactions_1_matches, reactions_2_matches = \
    get_subgraphs_overlap_max_results( graph_1, graph_2, results_subgraphs, results_matches, species_1, species_2, reactions_1, reactions_2)
    result = {}
    precision = 100. * len( matches_nodes_2) / float( len( graph_2.nodes()))
    recall  = 100. * len( matches_nodes_1) / float( len( graph_1.nodes()))
    result["node precision"] = precision
    result["node recall"] = recall
    if precision + recall == 0:    
        result["node f-score"] = 0
    else:
        result["node f-score"] = 2.0 * (precision * recall) / (precision + recall)
    
    precision = 100. * len( species_2_matches) / float( len( species_2))
    recall  = 100. * len( species_1_matches) / float( len( species_1))
    result["species precision"] = precision
    result["species recall"] = recall
    if precision + recall == 0:    
        result["species f-score"] = 0
    else:
        result["species f-score"] = 2.0 * (precision * recall) / (precision + recall)

    precision = 100. * len( reactions_2_matches) / float( len( reactions_2))
    recall  = 100. * len( reactions_1_matches) / float( len( reactions_1))
    result["reaction precision"] = precision
    result["reaction recall"] = recall
    if precision + recall == 0:    
        result["reaction f-score"] = 0
    else:
        result["reaction f-score"] = 2.0 * (precision * recall) / (precision + recall)

    precision = 100. * len( matches_edges_2) / float( len( graph_2.edges()))
    recall  =  100. * len( matches_edges_1) / float( len( graph_1.edges()))
    result["edge precision"] = precision
    result["edge recall"] = recall
    if precision + recall == 0:    
        result["edge f-score"] = 0
    else:
        result["edge f-score"] = 2.0 * (precision * recall) / (precision + recall)

    return result
    
def print_analysis_subgraphs_overlap_results( graph_1, graph_2, node_match, 
                                    matches_nodes_1, matches_nodes_2, matches_edges_1, matches_edges_2, 
                                    species_1_matches, species_2_matches, reactions_1_matches, reactions_2_matches,
                                    species_1 = None, species_2 = None,
                                    reactions_1 = None, reactions_2 = None):
    if not species_1:
        species_1 = set( filter_species( graph_1).nodes())
    if not species_2:
        species_2 = set( filter_species( graph_2).nodes())
    if not reactions_1:
        reactions_1 = set( filter_reactions( graph_1).nodes())
    if not reactions_2:
        reactions_2 = set( filter_reactions( graph_2).nodes())
    ## print results
    print( "{} {}/{}".format( node_match.__name__, graph_1.name, graph_2.name))
    precision = 100. * len( matches_nodes_2) / float( len( graph_2.nodes()))
    recall = 100. * len( matches_nodes_1) / float( len( graph_1.nodes()))
    f_score = 0.0
    if precision + recall > 0:
        f_score = 2. * precision * recall / (precision + recall)
    print( "%.2f & %.2f & %.2f node" % (precision, recall, f_score))

    precision = 100. * len( species_2_matches) / float( len( species_2))
    recall = 100. * len( species_1_matches) / float( len( species_1))
    f_score = 0.0
    if precision + recall > 0:
        f_score = 2. * precision * recall / (precision + recall)
    print( "%.2f & %.2f & %.2f species" % (precision, recall, f_score))

    precision = 100 * len( reactions_2_matches) / float( len( reactions_2))
    recall = 100 * len( reactions_1_matches) / float( len( reactions_1))
    f_score = 0.0
    if precision + recall > 0:
        f_score = 2. * precision * recall / (precision + recall)
    print( "%.2f & %.2f & %.2f reaction" % (precision, recall, f_score))

    precision = 100 * len( matches_edges_2) / float( len( graph_2.edges()))
    recall = 100 * len( matches_edges_1) / float( len( graph_1.edges()))
    f_score = 0.0
    if precision + recall > 0:
        f_score = 2. * precision * recall / (precision + recall)
    print( "%.2f & %.2f & %.2f edge" % (precision, recall, f_score))

    
def print_analysis_subgraphs_overlap_results_from_file( graph_1_name, 
                                                       graph_2_name, 
                                                       node_match,
                                                       edge_match = edge_match_exact,
                                                       prefix = "results/results-subgraphs-overlap-max"):
    # load file
    [graph_1, graph_2, results_subgraphs, results_matches] \
    = pickle.load( open( "%s__%s__%s__%s__%s.pickle" % ( prefix,graph_1_name, graph_2_name, node_match.__name__, edge_match.__name__), "rb"))
    # process results
    species_1 = set(filter_species( graph_1).nodes())
    species_2 = set(filter_species( graph_2).nodes())
    reactions_1 = set( filter_reactions( graph_1).nodes())
    reactions_2 = set( filter_reactions( graph_2).nodes())

    matches_nodes_1, matches_nodes_2, \
    matches_edges_1, matches_edges_2, \
    species_1_matches, species_2_matches, \
    reactions_1_matches, reactions_2_matches = \
    get_subgraphs_overlap_max_results( graph_1, graph_2, results_subgraphs, results_matches, \
                                   species_1 = species_1, species_2 = species_2, reactions_1 = reactions_1, reactions_2 = reactions_2)
    # print results
    print_analysis_subgraphs_overlap_results( graph_1, graph_2, node_match,
                                             matches_nodes_1, matches_nodes_2, matches_edges_1, matches_edges_2, 
                                             species_1_matches, species_2_matches, reactions_1_matches, reactions_2_matches,
                                             species_1, species_2, reactions_1, reactions_2)
                                   
def run_analysis_subgraphs_overlap( graph_1, 
                                    graph_2,
                                    node_match,
                                    edge_match = edge_match_exact,
                                    subgraphs_2 = None,
                                    species_1 = None, 
                                    species_2 = None, 
                                    reactions_1 = None, 
                                    reactions_2 = None,
                                    n_jobs = None,
                                    export_results = False,
                                    export_results_prefix = "results-subgraphs-overlap-max",
                                    file_names = None,
                                    print_results = True,
                                    ignore_existing = False):
    """ runs analysis for subgraphs """    
    print( "-----------")
    print( "%s: run_analysis_subgraphs_overlap %s/%s -- %s" % (now(), graph_1.name, graph_2.name, node_match.__name__))                                    
    if subgraphs_2 == None:
        subgraphs_2 = list( networkx.connected_component_subgraphs( graph_2))
    if not species_1:
        species_1 = set( filter_species( graph_1).nodes())
    if not species_2:
        species_2 = set( filter_species( graph_2).nodes())
    if not reactions_1:
        reactions_1 = set( filter_reactions( graph_1).nodes())
    if not reactions_2:
        reactions_2 = set( filter_reactions( graph_2).nodes())
        
    results_matches, results_subgraphs \
    = compute_subgraphs_overlap_max(  graph_1, graph_2, 
                                    node_match = node_match, 
                                    edge_match = edge_match_exact, 
                                    subgraphs_2 = subgraphs_2, 
                                    n_jobs = n_jobs, 
                                    export_results = export_results,
                                    export_results_prefix = export_results_prefix,
                                    file_names = file_names,
                                    ignore_existing = ignore_existing)
    if print_results:
        # process results
        matches_nodes_1, matches_nodes_2, \
        matches_edges_1, matches_edges_2, \
        species_1_matches, species_2_matches, \
        reactions_1_matches, reactions_2_matches = \
        get_subgraphs_overlap_max_results( graph_1, graph_2, results_subgraphs, results_matches, \
                                       species_1 = species_1, species_2 = species_2, reactions_1 = reactions_1, reactions_2 = reactions_2)
        # print results
        print_analysis_subgraphs_overlap_results( graph_1, graph_2, node_match,
                                                 matches_nodes_1, matches_nodes_2, matches_edges_1, matches_edges_2, 
                                                 species_1_matches, species_2_matches, reactions_1_matches, reactions_2_matches,
                                                 species_1, species_2, reactions_1, reactions_2)
        
    return results_matches, results_subgraphs

def run_analyses_subgraphs_overlap( graph_1, 
                                    graph_2, 
                                    node_match_fns,
                                    subgraphs_2 = None,
                                    export_results = False,
                                    export_results_prefix = "results-subgraphs-overlap-max",
                                    print_results = True):
    """ runs analysis for subgraphs  for multiple node_match_fns"""                                        
    print( "-----------")
    print( "run_analyses_subgraphs_overlap {}/{} -- {}".format( graph_1.name, graph_2.name, node_match_fns))
    
    if subgraphs_2 == None:
        subgraphs_2 = list( networkx.connected_component_subgraphs( graph_2))
    species_1 = set( filter_species( graph_1).nodes())
    species_2 = set( filter_species( graph_2).nodes())
    reactions_1 = set( filter_reactions( graph_1).nodes())
    reactions_2 = set( filter_reactions( graph_2).nodes())

    for node_match in node_match_fns:
        print("\n---")
        run_analysis_subgraphs_overlap( graph_1, graph_2, node_match, edge_match = edge_match_exact, 
                                    subgraphs_2 = subgraphs_2,
                                    species_1 = species_1, 
                                    species_2 = species_2, 
                                    reactions_1 = reactions_1, 
                                    reactions_2 = reactions_2,
                                    export_results = export_results,
                                    export_results_prefix = export_results_prefix)

########################################################################
########################################################################
## side-by-side graphviz

def _graphviz_label( n, graph, n_id_map = {}, id_prefix = "", participant = False, show_identifier = False):
    if graph.node[n].get("participants"):
        n_id_map[n] = id_prefix + n
        label = graph.node[n]["name"]
        if show_identifier:
            label += "\n" + n
        label = "<table>%s%s</table>" % ("<tr><td port=\"%s\"><b>%s</b></td></tr>" % (n, label),
                                   "".join([ _graphviz_label( p, graph, n_id_map, id_prefix = id_prefix, participant = True) for p in graph.node[n]["participant_ids"]]))
        if participant:
            return "<tr><td>%s</td></tr>" % label
        else:
            return "<%s>" % label
    elif graph.node[n]["type"] == "species":
        n_id_map[n] = id_prefix + n
        label = graph.node[n]["name"]
        if show_identifier:
            label += "\n" + n
        if participant:
            return "<tr><td port=\"%s\">%s</td></tr>" % (n, label)
        else:
            return label
    else:
        n_id_map[n] = n
        label = ", ".join( sbo_go_name(b) for b in graph.node[n]["bqbiol_is"])
        if show_identifier:
            label += "\n" + n
        return label
        
def _graphviz_add_node( n, graph, graphviz_graph, n_id_map = {}, label = None, show_identifier = False, **kwargs):
    """ adds a node top level (should not be participant of a complex) """
    if label == None and graph.node[n].get("participants"):
        label = _graphviz_label( n, graph, n_id_map, id_prefix = n + ":", show_identifier = show_identifier)
    else:
        label = _graphviz_label( n, graph, n_id_map, show_identifier = show_identifier)
    if graph.node[n].get("participants"): # has participants
        graphviz_graph.node( n, label = label, shape = "none", **kwargs)
    elif graph.node[n]["type"] == "species":
        graphviz_graph.node( n, label = label, shape = "rectangle", **kwargs)
    else:
        graphviz_graph.node( n, label = label, shape = "ellipse", **kwargs)

def _graphviz_add_edge( e, graph, graphviz_graph, n_id_map = {}, **kwargs):
    """ adds an edge to the graphviz graph """
    if (e[2] == "product" and not graph.node[e[0]]["type"] == "reaction") \
     or (e[2] != "product" and not graph.node[e[1]]["type"] == "reaction"):
         e = (e[1],e[0],e[2])
    e0 = e[0]
    e1 = e[1]
    if e0 in n_id_map:
        e0 = n_id_map[e0]
    if e1 in n_id_map:
        e1 = n_id_map[e1]
    if e[2] == "modifier":
        graphviz_graph.edge( e0, e1, arrowhead = "diamond", **kwargs)
    else:
        graphviz_graph.edge( e0, e1, **kwargs)        

def graphviz_graph( graph, file_name = "test.dot", view = True, show_identifier = False):
    """ renders a graph using dot"""
    import graphviz
    n_id_map = {}
    participant_complex_map = { id : n for n in graph.nodes() if graph.node[n].get("participant_ids") for id in graph.node[n].get("participant_ids")}
    top_nodes = set( graph.nodes()).difference( participant_complex_map.keys())
    graphviz_graph = graphviz.Digraph()
    [_graphviz_add_node( n, graph, graphviz_graph, n_id_map, show_identifier = show_identifier) for n in top_nodes];
    [_graphviz_add_edge( (e[0], e[1], e[2]["type"]), graph, graphviz_graph, n_id_map) for e in graph.edges( data = True)];
    
    graphviz_graph.render( file_name, view = view)

def get_top_complex( n, graph, participant_complex_map = {}):
    if participant_complex_map == {}:
        participant_complex_map = { id : n for n in graph.nodes() if graph.node[n].get("participant_ids") for id in graph.node[n].get("participant_ids")}
    if not (n in participant_complex_map):
        return n
    else:
        return get_top_complex( participant_complex_map[n], graph, participant_complex_map)
        
def graphviz_comparison_graph( graph_1, graph_2, m_nodes, m_edges, file_name = "test.dot", view = True, include_context_graph_1 = True, show_identifier = False):
    """ Creates a graph visualization visualizing a match (left, right) """
    import graphviz
    g = graphviz.Digraph()
    
    s1 = graphviz.Digraph( "cluster_1")
    s1.body.append( "\tlabel=\"%s\"" % graph_1.name)
    
    participant_complex_map_1 = { id : n for n in graph_1.nodes() if graph_1.node[n].get("participant_ids") for id in graph_1.node[n].get("participant_ids")}
    top_complexes_1 = [ get_top_complex( n, graph_1, participant_complex_map_1) for n in m_nodes.keys()]
    n_id_map_1 = {}
    
    [_graphviz_add_node( n, graph_1, s1, n_id_map_1, show_identifier = show_identifier) for n in top_complexes_1]
    [_graphviz_add_edge( e, graph_1, s1, n_id_map_1) for e in m_edges.keys()]
    if include_context_graph_1:
        context_edges = set([ sort_edge_signature((edge[0], edge[1], edge[2]["type"]), graph_1) for edge in graph_1.edges( m_nodes.keys(), data = True)]).difference( m_edges.keys())
        context_nodes = set( [ c for e in context_edges for c in e[:2]]).difference( top_complexes_1)
        # add nodes
        [_graphviz_add_node( get_top_complex( n, graph_1, participant_complex_map_1), graph_1, s1, n_id_map_1, color = "grey", show_identifier = show_identifier) for n in context_nodes]
        # add edges
        [_graphviz_add_edge( e, graph_1, s1, n_id_map_1, color = "grey") for e in context_edges]
        
    g.subgraph( s1)
    
    s2 = graphviz.Digraph( "cluster_2")
    s2.body.append( "\tlabel=\"%s\"" % graph_2.name)
    participant_complex_map_2 = { id : n for n in graph_2.nodes() if graph_2.node[n].get("participant_ids") for id in graph_2.node[n].get("participant_ids")}
    top_complexes_2 = [ get_top_complex( n, graph_2, participant_complex_map_2) for n in m_nodes.values()]
    n_id_map_2 = {}
    [_graphviz_add_node( n, graph_2, s2, n_id_map_2, show_identifier = show_identifier) for n in top_complexes_2]
    [_graphviz_add_edge( e, graph_2, s2, n_id_map_2) for e in m_edges.values()]
    g.subgraph( s2)
    
    for n1, n2 in m_nodes.iteritems():
        g.edge( n_id_map_1[n1], n_id_map_2[n2], dir = "none", style = "dotted", constraint = "false")
    g.render( file_name, view = view)
    
#graphviz_comparison_graph( graph_1, graph_2, m_nodes, m_edges)
#graphviz_comparison_graph( graph_1, graph_2, m_nodes, m_edges, include_context_graph_1 = False)

########################################################################
########################################################################
## graphviz color overlap

def _graphviz_label2( n, graph, n_id_map = {}, id_prefix = "",  matched_nodes = set(), participant = False, **kwargs):
    """ Generates colored labels condition on whether a node matched """
    # choose color
    color = kwargs["color"]
    fontcolor = kwargs["fontcolor"]
    fillcolor = kwargs["fillcolor"]
    show_identifier = kwargs["show_identifier"]
    if n in matched_nodes:
        color = kwargs["matched_color"]
        fontcolor = kwargs["matched_fontcolor"]
        fillcolor = kwargs["matched_fillcolor"]

    # handle different node types (with or without participants etc)
    if graph.node[n].get("participants"):
        n_id_map[n] = id_prefix + n
        label = "<table>%s%s</table>" \
            % ("<tr><td port=\"%s\" color=\"%s\" bgcolor=\"%s\"><font color=\"%s\"><b>%s</b></font></td></tr>" % (n, color, fillcolor, fontcolor, graph.node[n]["name"]),
               "".join([ _graphviz_label2( p, graph, n_id_map, id_prefix = id_prefix, participant = True, **kwargs) for p in graph.node[n]["participant_ids"]]))
        if participant:
            return "<tr><td>%s</td></tr>" % label
        else:
            return "<%s>" % label
    elif graph.node[n]["type"] == "species":
        n_id_map[n] = id_prefix + n
        if participant:
            return "<tr><td port=\"%s\" color=\"%s\" bgcolor=\"%s\"><font color=\"%s\">%s</font></td></tr>" % (n, color, fillcolor, fontcolor, graph.node[n]["name"])
        else:
            label = graph.node[n]["name"]
            if show_identifier:
                label += "\n" + n
            return label
    else:
        n_id_map[n] = n
        label = ", ".join( sbo_go_name(b) for b in graph.node[n]["bqbiol_is"])
        if show_identifier:
            label += "\n" + n
        return label
            
        
def _graphviz_add_node2( n, graph, graphviz_graph, n_id_map = {}, matched_nodes = set(),
                         **kwargs):

    # choose color
    color = kwargs["color"]
    fontcolor = kwargs["fontcolor"]
    fillcolor = kwargs["fillcolor"]
    if n in matched_nodes:
        color = kwargs["matched_color"]
        fontcolor = kwargs["matched_fontcolor"]
        fillcolor = kwargs["matched_fillcolor"]
    
    # compute label        
    if graph.node[n].get("participants"):
        label = _graphviz_label2( n, graph, n_id_map, id_prefix = n + ":", matched_nodes = matched_nodes, **kwargs)
    else:
        label = _graphviz_label2( n, graph, n_id_map, matched_nodes = matched_nodes, **kwargs)
    if graph.node[n].get("participants"): # has participants
        graphviz_graph.node( n, label = label, color = color, fontcolor = fontcolor, fillcolor = fillcolor, shape = "none")
    elif graph.node[n]["type"] == "species": # simple species
        graphviz_graph.node( n, label = label, shape = "rectangle", color = color, fontcolor = fontcolor, fillcolor = fillcolor, style = "filled")
    else: # simple reaction
        graphviz_graph.node( n, label = label, shape = "ellipse", color = color, fontcolor = fontcolor, fillcolor = fillcolor, style = "filled")
        
def graphviz_comparison_graph2( graph_1, 
                              matched_nodes = set(), 
                              matched_edges = set(), 
                              file_name = "test.dot", 
                              view = True, 
                              mode = "only_match", # can be only_match, context, all
                              matched_color = "red",
                              matched_fontcolor = "red",
                              matched_fillcolor = "white",
                              fontcolor = "grey",
                              color = "grey",
                              fillcolor = "white",
                              show_identifier = False):
    """ Visualization of matched nodes and edges using different color (single graph)"""
    import graphviz
    g = graphviz.Digraph()
        
    participant_complex_map_1 = { id : n for n in graph_1.nodes() if graph_1.node[n].get("participant_ids") for id in graph_1.node[n].get("participant_ids")}
    top_complexes_1 = [ get_top_complex( n, graph_1, participant_complex_map_1) for n in matched_nodes]
    n_id_map_1 = {}
    
    for n in top_complexes_1:
        _graphviz_add_node2( n, graph_1, g, n_id_map_1, matched_nodes, 
                            matched_color = matched_color, 
                            matched_fontcolor = matched_fontcolor, 
                            matched_fillcolor = matched_fillcolor, 
                            fontcolor = fontcolor, 
                            color = color, 
                            fillcolor = fillcolor, 
                            show_identifier = show_identifier)

    [_graphviz_add_edge( e, graph_1, g, n_id_map_1, color = matched_color) for e in matched_edges]

    if mode == "context":
        
        context_edges = set([ sort_edge_signature((edge[0], edge[1], edge[2]["type"]), graph_1) for edge in graph_1.edges( matched_nodes, data = True)]).difference( matched_edges)
        context_nodes_complexes = set([get_top_complex( n, graph_1, participant_complex_map_1) for n in set( [ c for e in context_edges for c in e[:2]]).difference( matched_nodes)]).difference(top_complexes_1)
        # add context nodes
        [_graphviz_add_node2( n, graph_1, g, n_id_map_1, color = color, fontcolor = fontcolor, fillcolor = fillcolor, show_identifier = show_identifier) for n in context_nodes_complexes]
        # add context edges
        [_graphviz_add_edge( e, graph_1, g, color = color) for e in set(context_edges).difference( matched_edges)]
    elif mode == "all":
        all_top_complexes = set( [ get_top_complex( n, graph_1, participant_complex_map_1) for n in set(graph_1.nodes()).difference( matched_nodes)])
        all_edges = set([ (edge[0], edge[1], edge[2]["type"]) for edge in graph_1.edges( data = True)]).difference( matched_edges)
        # add context nodes
        [_graphviz_add_node2( n, graph_1, g, n_id_map_1, color = color, fontcolor = fontcolor, fillcolor = fillcolor, show_identifier = show_identifier) for n in all_top_complexes]
        # add context edges
        [_graphviz_add_edge( e, graph_1, g, color = color) for e in all_edges]
        
        
    g.render( file_name, view = view)

#graphviz_comparison_graph2( graph_1, set(m_nodes.keys()), set(m_edges.keys()))
#graphviz_comparison_graph2( graph_1, set(m_nodes.keys()), set(m_edges.keys()), mode = "all")

########################################################################
########################################################################

def subgraph_overlap_graphviz( file_name = "TARGET__NLP-ANN__nm_name_clean_approx_OR_gene_id_intersect_AND_sbo_is_a__edge_match_exact--MAX.pickle"):
    """ Creates a single overlap graph from subgraph match results (color based) """
    import pickle
    
    [graph_1, graph_2, subgraphs_2, matches_list] = pickle.load( open( file_name, "rb"))
    edge_signatures_1 = edge_signatures( graph_1)
    edge_signatures_2 = edge_signatures( graph_2)
    all_matched_nodes_1 = set()
    all_matched_edges_1 = set()
    for subgraph_2, matches in zip( subgraphs_2, matches_list):
        for m in matches:
            all_matched_nodes_1.update( m.keys())
            reverse_m = { v: k for k, v in m.iteritems()}
            m_edges = subgraph_match_get_edges( subgraph_2, m , reverse_m, edge_signatures_1, edge_signatures_2)
            all_matched_edges_1.update( m_edges.keys())
    
    graphviz_comparison_graph2( graph_1, all_matched_nodes_1, all_matched_edges_1, mode = "all", file_name = file_name + ".dot")

def subgraph_overlaps_graphviz( input_file = "results/results-subgraphs-overlap-max__TARGET__NLP-ANN__nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_overlaps_sbo_is_a__edge_match_exact.pickle",
                              output_file_prefix = "results-subgraphs-overlap-max__TARGET__NLP-ANN__nm_name_clean_approx_OR_bqbiol_is_equal_AND_nm_bqbiol_is_overlaps_sbo_is_a__edge_match_exact",
                              include_context_graph_1 = True,
                              ignore_isolated_nodes = True,
                              graph_1 = None,
                              graph_2 = None,
                              show_identifier = False,
                              reactions_1 = None,
                              reactions_2 = None,
                              graph_1_reaction_txt_mapping = None,
                              graph_2_reaction_txt_mapping = None):
    """ Creates many overlap graph from subgraph match results (comparison graph left/right) """
    import pickle
    
    [graph_1_f, graph_2_f, subgraphs_2, matches_list] = pickle.load( open( input_file, "rb"))
    if graph_1 == None:
        graph_1 = graph_1_f
    if graph_2 == None:
        graph_2 = graph_2_f
    if reactions_1 == None:
        reactions_1 = set( filter_reactions( graph_1))
    if reactions_2 == None:
        reactions_2 = set( filter_reactions( graph_2))
    edge_signatures_1 = edge_signatures( graph_1)
    edge_signatures_2 = edge_signatures( graph_2)
    for i, subgraph_2, matches in zip( range(len(subgraphs_2)), subgraphs_2, matches_list):
        print( "Processing %i of %i" % (i, len(subgraphs_2)))
        if ignore_isolated_nodes and len(subgraph_2.nodes()) < 2:
            print( "Ignoring %i of %i" % (i, len(subgraphs_2)))
                
        else:
            for j, m_nodes in enumerate( matches):
                print( "Processing matches %i of %i" % (j, len(matches)))
                reverse_m_nodes = { v: k for k, v in m_nodes.iteritems()}
                m_edges = subgraph_match_get_edges( subgraph_2, m_nodes, reverse_m_nodes, edge_signatures_1, edge_signatures_2)
                output_file = "%s-%i-%i.dot" % ( output_file_prefix, i, j)
                print( "Exporting %s" % output_file)
                graphviz_comparison_graph( graph_1, 
                                          graph_2, 
                                          m_nodes, 
                                          m_edges, 
                                          file_name = output_file, 
                                          view = False, 
                                          show_identifier = show_identifier,
                                          include_context_graph_1 = include_context_graph_1)
                if graph_1_reaction_txt_mapping:
                    output_file = "%s-%i-%i-target.txt" % ( output_file_prefix, i, j)
                    print( "Exporting %s" % output_file)
                    m_reactions_1 = reactions_1.intersection( m_nodes.keys())
                    open( output_file, "wt").write( "\n".join( [graph_1_reaction_txt_mapping[r] for r in m_reactions_1 if r in graph_1_reaction_txt_mapping]))
                
                if graph_2_reaction_txt_mapping:
                    output_file = "%s-%i-%i-nlp.txt" % ( output_file_prefix, i, j)
                    print( "Exporting %s" % output_file)
                    m_reactions_2 = reactions_2.intersection( m_nodes.values())
                    open( output_file, "wt").write( "\n".join( [graph_2_reaction_txt_mapping[r] for r in m_reactions_2 if r in graph_2_reaction_txt_mapping]))
                                          
########################################################################
########################################################################
## overlap SBML
                                      
def _sbml_color_all( root, color_lines = "90000000", color_bounds = "00000000"):
    namespaces = {"cd" : "http://www.sbml.org/2001/ns/celldesigner", "sbml" : "http://www.sbml.org/sbml/level2"}
    for line in root.xpath("//cd:line", namespaces = namespaces):
        line.set( "color", color_lines)
    for paint in root.xpath("//cd:paint", namespaces = namespaces):
        paint.set( "color", color_bounds)

def _sbml_color_reaction( root, reaction_id, color = "ffff0000", width = "1.0"):
    """ colors the reaction links to reactant and product"""
    namespaces = {"cd" : "http://www.sbml.org/2001/ns/celldesigner", "sbml" : "http://www.sbml.org/sbml/level2"}
    lines = root.xpath("//sbml:reaction[@id='%s']/sbml:annotation/cd:line" % reaction_id, namespaces = namespaces)
    assert( len(lines) == 1)
    lines[0].set( "color", color)
    lines[0].set( "width", width)

def _sbml_color_reaction_modifier( root, reaction_id, modifier_id, color = "ffff0000", width = "1.0"):
    namespaces = {"cd" : "http://www.sbml.org/2001/ns/celldesigner", "sbml" : "http://www.sbml.org/sbml/level2"}
    lines = root.xpath("//sbml:reaction[@id='%s']/sbml:annotation/cd:listOfModification/cd:modification[@aliases='%s']/cd:line" % (reaction_id, modifier_id), namespaces = namespaces)
    if len(lines) == 1:
        lines[0].set( "color", color)
        lines[0].set( "width", width)
    else:
        print( "_sbml_color_reaction_modifier:Ignoring %s/%s" % (reaction_id, modifier_id))

def _sbml_color_species( root, species_id, color = "ffff0000"):
    namespaces = {"cd" : "http://www.sbml.org/2001/ns/celldesigner", "sbml" : "http://www.sbml.org/sbml/level2"}
    paints = root.xpath( "//cd:speciesAlias[@id='%s']//cd:paint" % species_id, namespaces = namespaces) \
        or root.xpath( "//cd:complexSpeciesAlias[@id='%s']//cd:paint" % species_id, namespaces = namespaces) 
    assert( len(paints) > 0)
    [ p.set( "color", color) for p in paints]


def subgraph_overlaps_sbml( graph_1, 
                        matches_nodes_1 = set(),
                        matches_edges_1 = set(), 
                        inn = 'mTORPathway-celldesigner.xml',
                        out = 'mTORPathway-celldesigner-color.xml',
                        background_color_bounds = "00000000",
                        background_color_lines = "000000",
                        matched_color = "FF00FF00",                        
                        matched_line_width = "2.0"):
    """ Visualization of matched species and reactions using different color (single graph)"""
    print( "sbml_color_matched:Loading %s" % inn)
    tree = lxml.etree.parse( inn);
    root = tree.getroot()

    print( "subgraph_overlaps_sbml:Coloring background")
    _sbml_color_all( root, color_bounds = background_color_bounds, color_lines = background_color_lines)
    # color species
    print( "subgraph_overlaps_sbml:Coloring matched species")
    for n in set( filter_species( graph_1).nodes()).intersection( matches_nodes_1):
        _sbml_color_species( root, n, color = matched_color)
    
    print( "subgraph_overlaps_sbml:Coloring matched reactions")
    matched_reactions = set( filter_reactions( graph_1).nodes()).intersection( matches_nodes_1)
    modifier_edges = filter( lambda e: e[2] == "modifier", matches_edges_1)
    matched_modifiers = { r : [e[0] for e in modifier_edges if e[1] == r] for r in matched_reactions}
    for r in matched_reactions:
        _sbml_color_reaction( root, r, color = matched_color, width = matched_line_width)
        for m in matched_modifiers[r]:
            _sbml_color_reaction_modifier( root, r, m, color = matched_color, width = matched_line_width)
    
    print( "subgraph_overlaps_sbml:Outputting %s" % out)
    tree.write( out, encoding='utf-8', xml_declaration = True) 

                                                                       
########################################################################
########################################################################
# initialize
def initialize():
    global SBO_NODES, GENE_MAP, SIMSTRING_DB
    print( "Initializing networkx_analysis.py")
    print( "Loading SBO")
    SBO_NODES = pickle.load( open( "sbo.pickle", "rb"))
    print( "Loading GENE_MAP")
    GENE_MAP = pickle.load( open( "gene_map.pickle", "rb"))
    print( "Loading SIMSTRING_DB")
    SIMSTRING_DB = simstring.reader( 'gene_list.simstring')
    SIMSTRING_DB.measure = simstring.cosine
    SIMSTRING_DB.threshold = 0.9

def load_pathway( name,
                 input_file,
                 output_file,
                 output_file_participant_graph = None,
                 output_file_w_participant_edges = None,
                 ending = ".xml",
                 pickle_graph = True,
                 prefix = ""):
    # load data
    print( "Loading %s" % (input_file))
    sbml = load_sbml( input_file)
    model = sbml.getModel();
    if model == None:
        print( "Error loading %s" % (input_file))
        return;
    
    graph, participant_graph, graph_w_participant_edges = create_graph( model, prefix = prefix)
    graph.name = name
    graph.source_file_name = input_file
    graph.file_name = output_file
    if pickle_graph == True:
        print( "Saving networkx as " + output_file)
        pickle_output_file = output_file
        
        networkx.write_gpickle( graph, pickle_output_file)
        if output_file_participant_graph:
            print( "Saving participant_graph networkx as " + output_file_participant_graph)
            networkx.write_gpickle( participant_graph, output_file_participant_graph)
        if output_file_w_participant_edges:
            print( "Saving graph with participant edges networkx as " + output_file_w_participant_edges)
            networkx.write_gpickle( graph_w_participant_edges, output_file_w_participant_edges)
  
    return graph, participant_graph, graph_w_participant_edges

########################################################################
########################################################################
# PROCESSING SIGNATURES

def run_analysis_bqbiol_is_signatures( bqbiol_is_1, bqbiol_is_2, 
                                      name_1 = "name_1", name_2 = "name_2", type = "species", 
                                      equal_fns = [operator.eq]):
    bqbiol_is_terms_set_1 = set( [b for t in bqbiol_is_1 for b in t])
    bqbiol_is_set_1 = set(bqbiol_is_1)
    bqbiol_is_terms_set_2 = set( [b for t in bqbiol_is_2 for b in t])
    bqbiol_is_set_2 = set(bqbiol_is_2)
    data = []
    
    res_1, res_2, precision, recall, f_score = analyse_set_overlap( bqbiol_is_terms_set_1, bqbiol_is_terms_set_2)
    print("%s:%s/%s:%s unique bqbiol_is terms equal: %.2f & %.2f & %.2f precision/recall/fscore" % (now(), name_1, name_2, type,
                                                                                                      precision, recall, f_score))
    data.append({ "graph_1" : name_1, "graph_2" : name_2, "unique" : True, "type" : type,
                  "reduction" : "bqbiol_is terms", "eq" : "eq",
                  "precision" : precision, "recall" : recall, "f-score" : f_score})
                  
    for eq_fun in equal_fns:
        res_1, res_2, precision, recall, f_score = analyse_set_overlap( bqbiol_is_set_1, bqbiol_is_set_2, eq_fun)
        print("%s:%s/%s:%s unique bqbiol_is signatures %s: %.2f & %.2f & %.2f precision/recall/fscore" % (now(), name_1, name_2, type, eq_fun.__name__,
                                                                                                          precision, recall, f_score))
        data.append({ "graph_1" : name_1, "graph_2" : name_2, "unique" : True, "type" : type,
                      "reduction" : "bqbiol_is signatures", "eq" : eq_fun.__name__,
                      "precision" : precision, "recall" : recall, "f-score" : f_score})
        
        res_1, res_2, precision, recall, f_score = analyse_list_overlap( bqbiol_is_1, bqbiol_is_2, eq_fun)
        print("%s:%s/%s:%s bqbiol_is signatures %s: %.2f & %.2f & %.2f precision/recall/fscore" % (now(), name_1, name_2, type, eq_fun.__name__,
                                                                                                   precision, recall, f_score))
        data.append({ "graph_1" : name_1, "graph_2" : name_2, "unique" : False, "type" : type,
                      "reduction" : "bqbiol_is signatures", "eq" : eq_fun.__name__,
                      "precision" : precision, "recall" : recall, "f-score" : f_score})
    return data
                  
def run_analysis_species_signatures( graph_1, graph_2, species_1 = None, species_2 = None):
    import pandas
    print("%s:%s/%s:run_analysis_species_signatures" % (now(), graph_1.name, graph_2.name))
    if species_1 == None:
        print("%s:%s/%s:run_analysis_species_signatures:filtering species graph_1" % (now(), graph_1.name, graph_2.name))
        species_1 = filter_species( graph_1)
    if species_2 == None:
        print("%s:%s/%s:run_analysis_species_signatures:filtering species graph_2" % (now(), graph_1.name, graph_2.name))
        species_2 = filter_species( graph_2)
    data = []
    print("%s:%s/%s:run_analysis_species_signatures:names" % (now(), graph_1.name, graph_2.name))
    for reduction_fun, equality_fn in zip( [clean_name, clean_name2, clean_name2], [operator.eq, operator.eq, name_approx_equal]):
    
        source_target = ([ reduction_fun( graph_1.node[n]["name"]) for n in species_1],
                         [ reduction_fun( graph_2.node[n]["name"]) for n in species_2])
        res_1, res_2, precision, recall, f_score = analyse_set_overlap( set(source_target[0]), set(source_target[1]), equality_fn)
        print("%s:%s/%s: species unique overlap %s/%s: %.2f & %.2f & %.2f precision/recall/fscore" % (now(), graph_1.name, graph_2.name, reduction_fun.__name__, equality_fn.__name__,
                                                    precision, recall, f_score))
        data.append({ "graph_1" : graph_1.name, "graph_2" : graph_2.name, "unique" : True, "type" : "species",
                      "reduction" : reduction_fun.__name__, "eq" : equality_fn.__name__,
                      "precision" : precision, "recall" : recall, "f-score" : f_score})
    
        res_1, res_2, precision, recall, f_score = analyse_list_overlap( source_target[0], source_target[1], equality_fn)
        print("%s:%s/%s: species overlap %s/%s: %.2f & %.2f & %.2f precision/recall/fscore" % (now(), graph_1.name, graph_2.name, reduction_fun.__name__, equality_fn.__name__,
                                                    precision, recall, f_score))
        data.append({ "graph_1" : graph_1.name, "graph_2" : graph_2.name, "unique" : False, "type" : "species", 
                      "reduction" : reduction_fun.__name__, "eq" : equality_fn.__name__,
                      "precision" : precision, "recall" : recall, "f-score" : f_score})
    # BQBIOL_IS
    print("%s:%s/%s:run_analysis_species_signatures:running bqbiol_is" % (now(), graph_1.name, graph_2.name))
    data.extend( run_analysis_bqbiol_is_signatures( bqbiol_is_1 = [ graph_1.node[n]["bqbiol_is"] for n in species_1], 
                                                   bqbiol_is_2 = [ graph_2.node[n]["bqbiol_is"] for n in species_2], 
                                                   name_1 = graph_1.name, name_2 = graph_2.name, type = "species", 
                                                   equal_fns = [ tuple_eq_empty_not_eq, tuple_overlaps]))
    data_p = pandas.DataFrame(data)
    return data_p

def run_analysis_reactions_signatures( graph_1, graph_2, reactions_1 = None, reactions_2 = None):
    import pandas
    print("%s:%s/%s:run_analysis_reactions_signatures" % (now(), graph_1.name, graph_2.name))
    if reactions_1 == None:
        print("%s:%s/%s:run_analysis_reactions_signatures:filtering reactions graph_1" % (now(), graph_1.name, graph_2.name))
        reactions_1 = filter_reactions( graph_1)
    if reactions_2 == None:
        print("%s:%s/%s:run_analysis_reactions_signatures:filtering reactions graph_2" % (now(), graph_1.name, graph_2.name))
        reactions_2 = filter_reactions( graph_2)
        
    reactions_w_bqbiol_is_1 = [ n for n in reactions_1 if graph_1.node[n]["bqbiol_is"]]
    reactions_w_bqbiol_is_2 = [ n for n in reactions_2 if graph_2.node[n]["bqbiol_is"]]
    bqbiol_is_1 = [ graph_1.node[n]["bqbiol_is"] for n in reactions_w_bqbiol_is_1]
    bqbiol_is_2 = [ graph_2.node[n]["bqbiol_is"] for n in reactions_w_bqbiol_is_2]
    data = run_analysis_bqbiol_is_signatures( bqbiol_is_1 = bqbiol_is_1,
                                                   bqbiol_is_2 = bqbiol_is_2, 
                                                   name_1 = graph_1.name, name_2 = graph_2.name, type = "reactions", 
                                                   equal_fns = [ tuple_eq_empty_not_eq, tuple_overlaps, tuple_overlaps_sbo_is_a])
    data_p = pandas.DataFrame(data)
    return data_p
    
def run_analysis_compartments_signatures( graph_1, graph_2, species_1 = None, species_2 = None):
    import pandas
    print("%s:%s/%s:run_analysis_compartments_signatures" % (now(), graph_1.name, graph_2.name))
    if species_1 == None:
        print("%s:%s/%s:run_analysis_compartments_signatures:filtering species graph_1" % (now(), graph_1.name, graph_2.name))
        species_1 = filter_species( graph_1)
    if species_2 == None:
        print("%s:%s/%s:run_analysis_compartments_signatures:filtering species graph_2" % (now(), graph_1.name, graph_2.name))
        species_2 = filter_species( graph_2)
    data = []
    print("%s:%s/%s:run_analysis_species_signatures names" % (now(), graph_1.name, graph_2.name))
    compartments_1 = set( [species_1.node[s]["compartment"] for s in species_1.nodes() if species_1.node[s]["compartment"]])
    compartments_2 = set( [species_2.node[s]["compartment"] for s in species_2.nodes() if species_2.node[s]["compartment"]])
    res_1, res_2, precision, recall, f_score = analyse_set_overlap( compartments_1, compartments_2)
    print("%s:%s/%s:compartment unique overlap eq: %.2f & %.2f & %.2f precision/recall/fscore" % (now(), graph_1.name, graph_2.name, precision, recall, f_score))
    data.append({ "graph_1" : graph_1.name, "graph_2" : graph_2.name, "unique" : True, "type" : "compartment",
                  "reduction" : "name", "eq" : "eq",
                  "precision" : precision, "recall" : recall, "f-score" : f_score})
    return pandas.DataFrame(data)


def run_analysis_signatures( graph_1, graph_2, export_file = None):
    import pandas
    print("%s:%s/%s:run_analysis_signatures:started" % (now(), graph_1.name, graph_2.name))
    data_s = run_analysis_species_signatures( graph_1, graph_2)
 
    print( "---")
    data_r = run_analysis_reactions_signatures( graph_1, graph_2)

    print( "---")
    data_c = run_analysis_compartments_signatures( graph_1, graph_2)
    data = data_s.append(data_r).append(data_c).reset_index()
    if export_file:
        print("%s:%s/%s:run_analysis_signatures:exporting %s" % (now(), graph_1.name, graph_2.name, export_file))
        data.to_pickle( export_file)
    return data
    

########################################################################
########################################################################
# PROCESSING NODE OVERLAP

def run_analyses_overlap( graph_1, graph_2,
                         species_nm = [ nm_name_equal, 
                                       nm_name_clean_equal, 
                                       nm_name_clean_approx, 
#                                       nm_gene_id_intersect, 
                                       nm_bqbiol_is_equal,
                                       nm_bqbiol_is_overlaps],
                         reactions_nm = [ nm_bqbiol_is_equal, 
                                         nm_bqbiol_is_overlaps, 
                                         nm_bqbiol_is_overlaps_sbo_is_a],
                         n_jobs = None,
                         export_results = True,
                         export_results_prefix = "results-nodes-overlap-max"):

    run_analyses_nodes_overlap_max( filter_species( graph_1),
                               filter_species( graph_2),
                               node_match_fns = species_nm,
                               n_jobs = n_jobs,
                               export_results = export_results,
                               export_results_prefix = export_results_prefix);

    #### REACTIONS OVERLAP
    run_analyses_nodes_overlap_max( filter_reactions( graph_1),
                               filter_reactions( graph_2),  
                               node_match_fns = reactions_nm,
                               n_jobs = n_jobs,
                               export_results = export_results,
                               export_results_prefix = export_results_prefix);
                               
                               