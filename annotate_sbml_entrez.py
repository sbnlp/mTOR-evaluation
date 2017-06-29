# -*- coding: utf-8 -*-
import libsbml
import os
import pickle
import urllib
import urllib2
import re
import sys

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

def query_gnat( text, url = 'http://textmining.ls.manchester.ac.uk:8081'):
    """ Query REST gnat API, returns hypotheses of Entrez gene ids ["string",[ids]]
        See http://textmining.ls.manchester.ac.uk:8081/?help for help"""
    # text ID (e.g., PubMed ID), text cross-reference (e.g., PubMed), entity type (e.g. gene, goterm), 
    # entity subtype (e.g. species of the gene, GO branch), entity candidate ID(s) [semi-colon separated] (e.g. 
    # Entrez gene ID, GO code), 0-based start position in the text, 0-based end position in 
    # the text, mention as found in the text, and a confidence score

    params = urllib.urlencode( { 'text': text } )
    response = urllib2.urlopen(url, params).read()
    if response_file_name:
        responses[text] = response
        pickle.dump( responses, open( response_file_name, "wb"))
    return response

def split_response( response):
    if response == "":
        return []
    else:
        response = response.strip( "\n").split("\n")
        response = [r.split("\t") for r in response]
        return [r[4].strip(";").split(";") for r in response]
    
def get_candidate_ids( text):
    """ Returns lists of lists of candidate gene ids """    
    return split_response( query_gnat( text))

def get_candidate_ids_flat( text):
    """ Returns lists of lists of candidate gene ids """    
    return set([id for id_list in split_response( query_gnat( text)) for id in id_list])
    

def check( value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1. """
    if value == None:
        err_msg = "Error processing libSBML returned a null value trying to {}".format( message)
        raise SystemExit( err_msg)
    elif type( value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error processing trying to {} libSBML returned error code {} : "{}"'.format( message, value, libsbml.OperationReturnValue_toString(value).strip())
            raise SystemExit(err_msg)
    else:
        return    

def find_or_add_cv_term( terms, type = libsbml.BIOLOGICAL_QUALIFIER, bio_type = libsbml.BQB_IS):
    term = None        
    if terms:
        for t in terms:
            if t.getQualifierType()  == type and t.getBiologicalQualifierType() == bio_type:
                term = t
                break;
    if term == None:
        term = libsbml.CVTerm();
        term.setQualifierType( type);
        term.setBiologicalQualifierType( bio_type);    
    return term

def clean_name2( name):
    return re.sub('[^a-zA-Z0-9-]', ' ', remove_prefixes( name.lower())).strip()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Add Entrez gene information (skip already present Entrez gene identifiers).')
    parser.add_argument( '--input',
                    action = "store", 
                    dest = "input",
                    default = None,
                    help = "CellDesigner SBML input file name.") 
    parser.add_argument( '--output',
                    action = "store",
                    dest = "output",
                    default = None,
                    help = "Pure SBML output file name.")
    cmd = parser.parse_args()
    if cmd.input is None:
        print( "Error. Please specify input.")
        sys.exit(-1)
    elif not os.path.isfile( cmd.input):
        print( "ERROR: input file %s does not exist.\n", cmd.input)
        sys.exit( 1)
    elif cmd.output is None:
        print( "Error. Please specify output.")
        sys.exit(-1)
    else:    

        print( "Loading sbml %s"% cmd.input)
        reader = libsbml.SBMLReader()
        document = reader.readSBML(  cmd.input)
        print( "Loaded %s (%s errors)".format(  cmd.input, document.getNumErrors()))
        model = document.getModel();
            
        # add bqb_is
        added = 0
        skipped = 0

        ## query GNAT for all entrez gene identifiers
        print( "Getting all species names")
        species_names = set( [clean_name2(species.getName()) for species in list(document.getModel().getListOfSpecies())])
        gnat_response = {}
        for n in species_names:
            print( "Retrieving Entrez Gene identifiers from GNAT for %s" % n)
            gnat_response[n] = get_candidate_ids(n)
        
        print( "Finished retrieval. Annotating ...")
        for idx, species in enumerate( list(document.getModel().getListOfSpecies())):
            sys.stdout.write( "\rProcessing {}/{}".format( idx + 1, document.getModel().getNumSpecies()))
            terms = species.getCVTerms()
            bqbiol_is_term = find_or_add_cv_term( terms)
            ids_list = responses[clean_name2(species.getName())]
            changed = False
            for id in [ id for ids in ids_list for id in ids]:
                uri = "urn:miriam:ncbigene:%08d" % int(id)
                if not uri in [bqbiol_is_term.getResourceURI(i) for i in range( bqbiol_is_term.getNumResources())]:
                    bqbiol_is_term.addResource( uri)
                    added += 1
                    changed = True
                else:
                    skipped += 1
            if changed:
                # just add - if it wasn't there yet it will be. if it was nothing changes
                species.addCVTerm( bqbiol_is_term)
                # check the namespaces
                assert( species.getAnnotation().getNumChildren() == 1) # only one top element
                rdf_annotation = species.getAnnotation().getChild(0)
                assert( rdf_annotation.getName() == "RDF") # top element is RDF
                rdf_annotation_namespaces = rdf_annotation.getNamespaces()
                rdf_annotation_namespaces_prefixes = [rdf_annotation_namespaces.getPrefix(i) for i in range(rdf_annotation_namespaces.getNumNamespaces())]
                prefixes = { "bqbiol" : "http://biomodels.net/biology-qualifiers/", 
                                "bqmodel" :"http://biomodels.net/model-qualifiers/"}
                for prefix in set(prefixes.keys()).difference(rdf_annotation_namespaces_prefixes):
                    rdf_annotation_namespaces.add( prefixes[prefix], prefix)
                    
        print( "\nFinished processing sbml. Added {}, already present {} entrez gene ids.".format( added, skipped))
            
        output_file_name = cmd_args.output
        if output_file_name is None:
            output_file_name = cmd_args.input_file_name
        print "Outputting {}".format( output_file_name)
        libsbml.writeSBMLToFile( document, output_file_name)