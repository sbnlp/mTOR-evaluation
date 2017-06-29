#!/usr/bin/python
# -*- coding: utf-8 -*-
from xml.dom import minidom
import cgi
import libsbml
import sys

CELLDESIGNER_TYPE_REFERENCE = {
    'GENE': 'geneReference',
    'RNA': 'rnaReference',
    'PROTEIN': 'proteinReference',
    'ANTISENSE_RNA': 'antisensernaReference'}

STANDOFF_ENTITY_TO_SBO_MAPPING = {  # informational molecule segment
                                    # non-covalent complex
                                    # polypeptide chain
                                    # deoxyribonucleic acid
                                    # deoxyribonucleic acid
                                    # ribonucleic acid
                                    # ribonucleic acid
                                    # simple chemical
                                    # simple chemical
                                    # simple chemical
                                    # simple chemical
    'gene': 'SBO:0000354',
    'complex': 'SBO:0000253',
    'protein': 'SBO:0000252',
    'dna': 'SBO:0000251',
    'dnaregion': 'SBO:0000251',
    'rna': 'SBO:0000250',
    'rnaregion': 'SBO:0000250',
    'smallmolecule': 'SBO:0000247',
    'simple_molecule': 'SBO:0000247',
    'ion': 'SBO:0000247',
    'drug': 'SBO:0000247',
    'phenotype': 'SBO:0000358'}

STANDOFF_EVENT_TO_SBO_MAPPING = {  # degradation
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
    'ubiquitination': 'SBO:0000224',
    'STATE_TRANSITION': 'SBO:0000182',
    'KNOWN_TRANSITION_OMITTED': 'SBO:0000182',
    'UNKNOWN_TRANSITION': 'SBO:0000182',
    'CATALYSIS': 'SBO:0000172',
    'UNKNOWN_CATALYSIS': 'SBO:0000172',
    'INHIBITION': 'SBO:0000169',
    'UNKNOWN_INHIBITION': 'SBO:0000169',
    'TRANSPORT': 'SBO:0000185',
    'HETERODIMER_ASSOCIATION': 'SBO:0000297',
    'DISSOCIATION': 'SBO:0000180',
    'TRUNCATION': 'SBO:0000180',
    'TRANSCRIPTIONAL_ACTIVATION': 'SBO:0000170',
    'TRANSCRIPTIONAL_INHIBITION': 'SBO:0000169',
    'TRANSLATIONAL_ACTIVATION': 'SBO:0000170',
    'TRANSLATIONAL_INHIBITION': 'SBO:0000169',
    'TRANSCRIPTION': 'SBO:0000183',
    'TRANSLATION': 'SBO:0000184'}

STANDOFF_EVENT_TO_SBO_MAPPING = dict((k.lower(), v) for (k, v) in STANDOFF_EVENT_TO_SBO_MAPPING.iteritems())

# mapping of general reaction components to sbo term

GENERIC_REACTION_SBO_MAPPING = {
    'reactant': 'SBO:0000010',
    'product': 'SBO:0000011',
    'modifier': 'SBO:0000019',
    'activator': 'SBO:0000021',
    'inhibitor': 'SBO:0000020'}

def add_cvterm(term):
    controlled_vocab = libsbml.CVTerm()
    controlled_vocab.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)
    controlled_vocab.setBiologicalQualifierType(libsbml.BQB_IS)
    controlled_vocab.addResource(term)
    return controlled_vocab

def add_note(note, species):
    """ Adds a note to species (wraps the note in <p></p> and escapes the text) """
    species.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(note)))
    return species

def add_annotation_complex(model, model_id, participants):

    species = model.getSpecies(model_id)
    xmlns = libsbml.XMLNamespaces()
    xmlns.add('http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'rdf')
    rdf_triple = libsbml.XMLTriple( 'RDF',
                                   'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 
                                   'rdf')
    rdf_token = libsbml.XMLToken( rdf_triple, 
                                 libsbml.XMLAttributes(),
                                 xmlns)
    annotation = libsbml.XMLNode(rdf_token)
    if species:
        participants_xml_triple = libsbml.XMLTriple( 'Participants',
                                                    'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 
                                                    'rdf')
        participants_xml_token = \
            libsbml.XMLToken(participants_xml_triple,
                             libsbml.XMLAttributes())
        participants_xml_node = libsbml.XMLNode(participants_xml_token)
        participant_xml_triple = libsbml.XMLTriple('Participant',
                                                   'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 
                                                   'rdf')
        for iii in participants:
            resource_att = libsbml.XMLAttributes()
            resource_att.add('participant', str(iii))
            subject_token = libsbml.XMLToken(participant_xml_triple,
                    resource_att)
            subject_token.setEnd()
            participants_xml_node.addChild(libsbml.XMLNode(subject_token))
        annotation.addChild(participants_xml_node)
    species.appendAnnotation(annotation)

def correct_species_name(species):
    if bool(species['modifications']):
        new_name = species['name']
        for modification in species['modifications']:
            new_name = species['modifications'][modification] + ' ' \
                + new_name
        species['newname'] = new_name
    return species

def get_complex_to_species_link( species, celldesigner_complex_alias):
    for i in celldesigner_complex_alias:
        if celldesigner_complex_alias[i]['species'] == species:
            return i

def correct_species_alias( species_alias, speciesIdentity, cell_designer_species_alias):
    if species_alias in cell_designer_species_alias.keys() \
        and cell_designer_species_alias[species_alias]['species'] \
        in speciesIdentity:
        if cell_designer_species_alias[species_alias]['activity'] \
            == 'active':
            if 'newname' \
                in correct_species_name(speciesIdentity[cell_designer_species_alias[species_alias]['species']]).keys():
                return 'activated ' \
                    + correct_species_name(speciesIdentity[cell_designer_species_alias[species_alias]['species']])['newname']
            else:
                return 'activated ' \
                    + correct_species_name(speciesIdentity[cell_designer_species_alias[species_alias]['species']])['name']
        else:
            if 'newname' \
                in correct_species_name(speciesIdentity[cell_designer_species_alias[species_alias]['species']]).keys():
                return correct_species_name(speciesIdentity[cell_designer_species_alias[species_alias]['species']])['newname']
            else:
                return correct_species_name(speciesIdentity[cell_designer_species_alias[species_alias]['species']])['name']

def add_compartment( compartment_id, model, compartmentIdentity):
    compartment_ref = model.getCompartment(str(compartment_id))
    if compartment_ref == None:
        compartment_ref = model.createCompartment()
        compartment_ref.setId(str(compartment_id))
        compartment_ref.setName(str(compartmentIdentity[compartment_id]['name']))
        compartment_ref.setConstant(True)
        compartment_ref.setSize(1)
        compartment_ref.setUnits('volume')
    return compartment_ref


def export_pure_sbml( input_file = 'mTORPathway-celldesigner.xml', output_file = 'mTORPathway-celldesigner_out.xml'):

    with open( input_file) as cellDesignerFile:
        cd_xml = cellDesignerFile.read()
    
    cd_xml_parsed = minidom.parseString(cd_xml)
    
    reactions = cd_xml_parsed.getElementsByTagName('reaction')
    species = cd_xml_parsed.getElementsByTagName('species')
    
    cell_designer_species = \
        cd_xml_parsed.getElementsByTagName('celldesigner:listOfIncludedSpecies')[0].getElementsByTagName('celldesigner:species')
    speciesAlias = \
        cd_xml_parsed.getElementsByTagName('celldesigner:speciesAlias')
    compartments = cd_xml_parsed.getElementsByTagName('listOfCompartments')[0].getElementsByTagName('compartment')
    complexAlias = \
        cd_xml_parsed.getElementsByTagName('celldesigner:complexSpeciesAlias')
    
    compartmentIdentity = {}
    for compartment in compartments:
        compartmentIdentity[compartment.attributes['id'].value] = {}
        if compartment.hasAttribute('name'):
            compartmentIdentity[compartment.attributes['id'].value]['name'] = compartment.attributes['name'].value
        else:
            compartmentIdentity[compartment.attributes['id'].value]['name'] = compartment.attributes['id'].value
    
    speciesIdentity = {}
    for specie in species:
        speciesIdentity[specie.attributes['id'].value] = {}
    
        # speciesIdentity[specie.attributes["id"].value]["metaid"] = specie.attributes["metaid"].value
    
        speciesIdentity[specie.attributes['id'].value]['name'] = \
            specie.attributes['name'].value
        speciesIdentity[specie.attributes['id'].value]['compartment'] = \
            specie.attributes['compartment'].value
        speciesIdentity[specie.attributes['id'].value]['Class'] = \
            specie.getElementsByTagName('celldesigner:speciesIdentity'
                                        )[0].getElementsByTagName('celldesigner:class'
                )[0].childNodes[0].data
        if speciesIdentity[specie.attributes['id'].value]['Class'] \
            in CELLDESIGNER_TYPE_REFERENCE.keys():
            speciesIdentity[specie.attributes['id'].value]['reference'] = \
                specie.getElementsByTagName('celldesigner:speciesIdentity'
                    )[0].getElementsByTagName('celldesigner:'
                    + CELLDESIGNER_TYPE_REFERENCE[speciesIdentity[specie.attributes['id'
                    ].value]['Class']])[0].childNodes[0].data
        speciesIdentity[specie.attributes['id'].value]['modifications'] = {}
        if specie.getElementsByTagName('celldesigner:listOfModifications'):
            if specie.getElementsByTagName('celldesigner:listOfModifications'
                    )[0].getElementsByTagName('celldesigner:modification'):
                modifications = \
                    specie.getElementsByTagName('celldesigner:listOfModifications'
                        )[0].getElementsByTagName('celldesigner:modification'
                        )
                for modification in modifications:
                    speciesIdentity[specie.attributes['id'
                                    ].value]['modifications'
                            ][modification.attributes['residue'].value] = \
                        modification.attributes['state'].value
    
    for specie in cell_designer_species:
        if specie.attributes['id'].value not in speciesIdentity:
            speciesIdentity[specie.attributes['id'].value] = {}
            speciesIdentity[specie.attributes['id'].value]['name'] = \
                specie.attributes['name'].value
            speciesIdentity[specie.attributes['id'].value]['Class'] = \
                specie.getElementsByTagName('celldesigner:speciesIdentity')[0].getElementsByTagName('celldesigner:class')[0].childNodes[0].data
            if speciesIdentity[specie.attributes['id'].value]['Class'] \
                in CELLDESIGNER_TYPE_REFERENCE.keys():
                speciesIdentity[specie.attributes['id'].value]['reference'] = \
                    specie.getElementsByTagName('celldesigner:speciesIdentity')[0].getElementsByTagName('celldesigner:' + CELLDESIGNER_TYPE_REFERENCE[speciesIdentity[specie.attributes['id'].value]['Class']])[0].childNodes[0].data
            speciesIdentity[specie.attributes['id'].value]['modifications'] = {}
            if specie.getElementsByTagName('celldesigner:listOfModifications'):
                if specie.getElementsByTagName('celldesigner:listOfModifications')[0].getElementsByTagName('celldesigner:modification'):
                    modifications = specie.getElementsByTagName('celldesigner:listOfModifications')[0].getElementsByTagName('celldesigner:modification')
                    for modification in modifications:
                        speciesIdentity[specie.attributes['id'
                                        ].value]['modifications'
                                ][modification.attributes['residue'].value] = \
                            modification.attributes['state'].value
    
    
    # ....else:
    # ........print specie.attributes["id"].value
    
    cell_designer_species_alias = {}
    for specie in speciesAlias:
        cell_designer_species_alias[specie.attributes['id'].value] = {}
        cell_designer_species_alias[specie.attributes['id'].value]['species'] = specie.attributes['species'].value
        cell_designer_species_alias[specie.attributes['id'].value]['activity'] = \
            specie.getElementsByTagName( 'celldesigner:activity')[0].childNodes[0].data
    
    celldesigner_complex_alias = {}
    for specie in complexAlias:
        celldesigner_complex_alias[specie.attributes['id'].value] = {}
        celldesigner_complex_alias[specie.attributes['id'].value]['players'] = []
        celldesigner_complex_alias[specie.attributes['id'].value]['species'] = specie.attributes['species'].value
        celldesigner_complex_alias[specie.attributes['id'].value]['activity'] = specie.getElementsByTagName('celldesigner:activity')[0].childNodes[0].data
    
    for specie in speciesAlias:
        if 'complexSpeciesAlias' in [item[0] for item in specie.attributes.items()]:
            celldesigner_complex_alias[specie.attributes['complexSpeciesAlias'].value]['players'].append(specie.attributes['id'].value)
    
    for specie in complexAlias:
        if 'complexSpeciesAlias' in [item[0] for item in specie.attributes.items()]:
            celldesigner_complex_alias[specie.attributes['complexSpeciesAlias'].value]['players'].append(specie.attributes['id'].value)
    
    reactionIdentity = {}
    for reaction in reactions:
        reactionIdentity[reaction.attributes['id'].value] = {}
    
        # reactionIdentity[reaction.attributes["id"].value]["metaid"] = reaction.attributes["metaid"].value
    
        reactionIdentity[reaction.attributes['id'].value]['reversible'] = \
            reaction.attributes['reversible'].value
        reactionIdentity[reaction.attributes['id'].value]['reactionType'] = \
            reaction.getElementsByTagName('celldesigner:reactionType')[0].childNodes[0].data
        reactants = {}
        products = {}
        modifiers = {}
        if reaction.getElementsByTagName('listOfReactants'):
            for reactant in reaction.getElementsByTagName('listOfReactants')[0].getElementsByTagName('speciesReference'):
                reactants[reactant.attributes['species'].value] = ''
        if reaction.getElementsByTagName('listOfProducts'):
            for product in reaction.getElementsByTagName('listOfProducts')[0].getElementsByTagName('speciesReference'):
                products[product.attributes['species'].value] = ''
        if reaction.getElementsByTagName('listOfModifiers'):
            for modifier in reaction.getElementsByTagName('listOfModifiers')[0].getElementsByTagName('modifierSpeciesReference'):
                modifiers[modifier.attributes['species'].value] = ''
        if reaction.getElementsByTagName('celldesigner:listOfModification'):
            listOfModifications = \
                reaction.getElementsByTagName('celldesigner:listOfModification')[0].getElementsByTagName('celldesigner:modification')
            for modification in listOfModifications:
                modifiers[modification.attributes['modifiers'].value] = \
                    [modification.attributes['type'].value,
                     modification.attributes['aliases'].value]
        if reaction.getElementsByTagName('celldesigner:baseReactants'):
            listOfModifications = \
                reaction.getElementsByTagName('celldesigner:baseReactants')[0].getElementsByTagName('celldesigner:baseReactant')
            for modification in listOfModifications:
                reactants[modification.attributes['species'].value] = \
                    modification.attributes['alias'].value
        if reaction.getElementsByTagName('celldesigner:baseProducts'):
            listOfModifications = \
                reaction.getElementsByTagName('celldesigner:baseProducts')[0].getElementsByTagName('celldesigner:baseProduct')
            for modification in listOfModifications:
                products[modification.attributes['species'].value] = \
                    modification.attributes['alias'].value
        reactionIdentity[reaction.attributes['id'].value]['reactants'] = \
            reactants
        reactionIdentity[reaction.attributes['id'].value]['products'] = \
            products
        reactionIdentity[reaction.attributes['id'].value]['modifiers'] = \
            modifiers
    
    document = libsbml.SBMLDocument(2, 4)
    model = document.createModel()
    
    for i in compartmentIdentity:
        add_compartment(i, model, compartmentIdentity)
    
    listofcell_designer_species = []
    
    for i in cell_designer_species_alias:
        if cell_designer_species_alias[i]['species'] in speciesIdentity \
            and speciesIdentity[cell_designer_species_alias[i]['species']]['Class'].lower() \
            in STANDOFF_ENTITY_TO_SBO_MAPPING:
            species = model.createSpecies()
            species.setId(str(i))
            species.setMetaId('metaid_0000' + str(i))
            listofcell_designer_species.append(str(i))
            species.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(cell_designer_species_alias[i]['activity'])))
            for j in \
                speciesIdentity[cell_designer_species_alias[i]['species']]['modifications']:
                species.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(speciesIdentity[cell_designer_species_alias[i]['species']]['modifications'][j])))
    
            species.setName( str( correct_species_alias( i, speciesIdentity, cell_designer_species_alias)))
    
            species.setSBOTerm(str(STANDOFF_ENTITY_TO_SBO_MAPPING[speciesIdentity[cell_designer_species_alias[i]['species']]['Class'].lower()]))
            if 'compartment' \
                in speciesIdentity[cell_designer_species_alias[i]['species']]:
                species.setCompartment(str(speciesIdentity[cell_designer_species_alias[i]['species']]['compartment']))
            else:
                species.setCompartment('default')

    for i in celldesigner_complex_alias:
        if celldesigner_complex_alias[i]['species'] in speciesIdentity \
            and speciesIdentity[celldesigner_complex_alias[i]['species']]['Class'].lower() \
            in STANDOFF_ENTITY_TO_SBO_MAPPING:
            species = model.createSpecies()
            species.setId(str(i))
            species.setMetaId('metaid_0000' + str(i))
            listofcell_designer_species.append(str(i))
            add_annotation_complex( model, str(i), celldesigner_complex_alias[i]['players'])
            for j in speciesIdentity[celldesigner_complex_alias[i]['species']]['modifications']:
                species.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(speciesIdentity[celldesigner_complex_alias[i]['species']]['modifications'][j])))
            species.setName(str(speciesIdentity[celldesigner_complex_alias[i]['species']]['name']))
            if 'compartment' \
                in speciesIdentity[celldesigner_complex_alias[i]['species']]:
                species.setCompartment(str(speciesIdentity[celldesigner_complex_alias[i]['species']]['compartment']))
            else:
                species.setCompartment('default')
            species.setSBOTerm(str(STANDOFF_ENTITY_TO_SBO_MAPPING[speciesIdentity[celldesigner_complex_alias[i]['species']]['Class'].lower()]))
    
    for i in reactionIdentity:
        reaction = model.createReaction()
        reaction.setName(str(i))
        reaction.setId(str(i))
        if reactionIdentity[i]['reversible'].upper() == 'TRUE':
            reaction.setReversible(True)
        else:
            reaction.setReversible(False)
        if STANDOFF_EVENT_TO_SBO_MAPPING.get(reactionIdentity[i]['reactionType'].lower()) \
            and (STANDOFF_EVENT_TO_SBO_MAPPING[reactionIdentity[i]['reactionType'].lower()])[0:3] == 'SBO':
            reaction.setSBOTerm(STANDOFF_EVENT_TO_SBO_MAPPING[reactionIdentity[i]['reactionType'].lower()])
        reaction.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(str(reactionIdentity[i]['reactionType']))))
        for j in reactionIdentity[i]['reactants']:
            if model.getSpecies(str(reactionIdentity[i]['reactants'][j])):
                reactant_ref = reaction.createReactant()
                reactant_ref.setSpecies(str(reactionIdentity[i]['reactants'][j]))
            elif model.getSpecies( str( get_complex_to_species_link( str(j), celldesigner_complex_alias))):
                reactant_ref = reaction.createReactant()
                reactant_ref.setSpecies(str(get_complex_to_species_link(str(j), celldesigner_complex_alias)))
        for j in reactionIdentity[i]['products']:
            if str(reactionIdentity[i]['products'][j]) \
                in listofcell_designer_species:
                product_ref = reaction.createProduct()
                product_ref.setSpecies(str(reactionIdentity[i]['products'][j]))
            elif str(get_complex_to_species_link( str(j), celldesigner_complex_alias)) \
                in listofcell_designer_species:
                product_ref = reaction.createProduct()
                product_ref.setSpecies(str(get_complex_to_species_link(str(j),celldesigner_complex_alias)))
        for j in reactionIdentity[i]['modifiers']:
            if reactionIdentity[i]['modifiers'][j]:
                if model.getSpecies(str(list(reactionIdentity[i]['modifiers'][j])[1])):
                    modifier_ref = reaction.createModifier()
                    modifier_ref.setSpecies(str(list(reactionIdentity[i]['modifiers'][j])[1]))
                    modifier_ref.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(str(list(reactionIdentity[i]['modifiers'][j])[0]))))
                    modifier_ref.appendNotes('<p xmlns="http://www.w3.org/1999/xhtml">{0}</p>'.format(cgi.escape(str(STANDOFF_EVENT_TO_SBO_MAPPING[reactionIdentity[i]['modifiers'][j][0].lower()]))))
            elif model.getSpecies(str(get_complex_to_species_link(str(j),celldesigner_complex_alias))):
                modifier_ref = reaction.createModifier()
                modifier_ref.setSpecies(str(get_complex_to_species_link(str(j),celldesigner_complex_alias)))
    
    libsbml.writeSBMLToFile(document, output_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
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
        export_pure_sbml( cmd.input, cmd.output)
