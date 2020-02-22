import time
import xml.etree.cElementTree as ET
import requests
from bs4 import BeautifulSoup
import libchebipy
import csv
namespace = '{http://uniprot.org/uniprot}'
## desired types = name, fullName, gene name, organism name, dbRef NCBI ID,
desired_types = {0: "name", 1: "fullName", 2: "gene_name", 3: "organism_name", 4: "NCBI", 5: "taxon", 6:"sequence", 7: "ChEBI", 8: "keyword"}
desired_tags = ["name","fullName", "dbReference", "taxon", "sequence", "reaction", "keyword"]
desired_elements = [] ### this is a list of all the elements, each element is a temp_entry which is also a list.

def parser(file_path):      
    path = []
    temp_elem = ["", "", "", "", "", "", "", "", "", None, ""]
    temp_chebis = []
    temp_chebi = []
    temp_taxon =[]
    temp_keywords = []
    for event, elem in ET.iterparse(file_path, events=("start", "end")):
        #print(event)
        tag = elem.tag.strip().replace(namespace, "")#.encode("utf-8")
        #print(tag)
        if event == 'start': ## if there is a new tag  
            ##add it to the path history
            path.append(tag)
            element_type = ""
            if tag in desired_tags:
                print(tag)
                print(elem.text)
                if  tag == "name" and "gene" in path:
                     element_type = "gene_" + tag
                elif tag == "name" and "protein" in path:
                     element_type = "protein_" + tag
                elif tag == "name" and "organism" in path and elem.get('type') == "scientific":
                     element_type = "organism_" + tag  
                elif tag == "dbReference" and "reaction" in path:
                    if elem.get('type') == "ChEBI":
                        temp_chebi.append(elem.get('id'))
                        element_type = "ChEBI"
                elif tag == "dbReference":
                    if elem.get('type') == "NCBI Taxonomy": 
                        element_type = "NCBI"
                elif tag == "taxon":
                    temp_taxon.append(elem.text)
                    element_type = "taxon"
                elif tag == "keyword":
                    temp_keywords.append(elem.text)
                    element_type = "keyword"
                elif tag == "name" and not elem.attrib:
                     element_type = "name"
                elif tag == "fullName":
                     element_type = tag  
                elif tag == "sequence":
                     element_type = "sequence"
            ### process further if the element_typeis actually wanted
            if element_type in desired_types.values():
                if element_type != "ChEBI" and element_type!= "taxon" and element_type !="keyword":   
                    index = list(desired_types.values()).index(element_type)
                    ## get the element_typeid if its dbReference
                    if  element_type == "NCBI":
                        temp_elem[index] = elem.get('id')
                    ### else if the content is contained in its text
                    elif elem.text:
                        temp_elem[index] = elem.text.strip()
                        #print(elem.text.strip())
                    #print(element, ' equal to ', temp_elem[index])
            elem.clear()
        elif event == 'end': ### end tag signals that a tag's contents r ready to be processed 
            #end of a singular chebi pairing
            if tag == "reaction":
                if len(temp_chebi) > 0:
                    temp_chebi = update_chebi(temp_chebi)
                    #print('reaction with chebis: ', temp_chebi)
                    temp_chebis.append(temp_chebi)
                    #print('reactions : ' , reactions)
                    temp_chebi = []   
            # end of all chebi streams, taxon, and keywords
            if tag == "entry":
                temp_elem[5] = temp_taxon
                #print('taxonomy : ', temp_taxon)
                temp_taxon = []
                #if there were one or more reactions
                if len(temp_chebis) > 0:
                    temp_elem[7] = temp_chebis
                    #clear all reactions
                    temp_chebis = []   
                if len(temp_keywords) > 0:
                    temp_elem[8] = temp_keywords
                    #clear all kwds
                    temp_keywords = []  
                desired_elements.append(temp_elem)
                #reset entry to null
                temp_elem = ["", "", "", "", "", "", "", "", "", None, ""]            
            path.pop()
            elem.clear()
def update_chebi(chebi_list):
    smiles = []
    for code in chebi_list: ### request smile and add to smiles[]
        smile = libchebipy.ChebiEntity(code).get_smiles()
        smiles.append(smile)
    return smiles

start = time.time()
parser("uniprot_trembl.xml")

#noe = 0
# for element in desired_elements:
#     noe +=1
#     print(' element # ', noe, ' name : ',  element[0], ' taxon :', element[6])
with open('2020-idk.csv', 'w+', newline = '') as file:
    writer = csv.writer(file)
    writer.writerow(desired_types.values())
    for element in desired_elements:
        print(element)
        writer.writerow(element)    
print('Time elapsed: ', (time.time() - start) / 60, ' minutes')    