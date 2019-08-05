# run by typing 
# python downloader.py
#
# Willy Bruhn, 24.6.19
#
# automatically downloads from
# https://www.drugbank.ca/
# redirected to here 
# https://www.rcsb.org/
# the pdb-files
#
#-----------------------------------------------------------

import urllib2
import re
import os  
from shutil import copyfile
from distutils.dir_util import copy_tree

from os import listdir
from os.path import isfile, join
import json



def create_dir(path):
    if not os.path.isdir(path):
        try:  
            os.mkdir(path)
        except OSError:  
            print ("Creation of the directory %s failed" % path)
        else:  
            print ("Successfully created the directory %s " % path)


def get_files_in_dir(path):
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    return(onlyfiles)


def get_all_drugs_for_target(target_name):
    polypeptides_url = 'https://www.drugbank.ca/polypeptides/'

    target_url = polypeptides_url + target_name
    target_url_file = target_name + '.html'



    if os.path.isfile(target_url_file):
        print(target_url + " exists")
    else:
        print("fetching data from " + target_url + " ...")
        response = urllib2.urlopen(target_url)
        webContent = response.read()
        f = open(target_url_file, "w")
        f.write(webContent)
        f.close()
        

    webContent = open(target_url_file, "r")
    webContentString = webContent.read()


    gene_name_search = re.search(r"Gene Name<\/dt><dd class=\"col-md-10 col-sm-8\">(.*)<\/dd><dt class=\"col-md-2 col-sm-4\">Organism",webContentString)

    gene_name = gene_name_search.group(1)
    print(gene_name)


    Organism_search = re.search(r"Organism<\/dt><dd class=\"col-md-10 col-sm-8\">(.*)<\/dd><dt class=\"col-md-2 col-sm-4\">Amino acid sequence",webContentString)

    organism_name = Organism_search.group(1)
    print(organism_name)

    # Now get the Drugs fitting to this target
    Drug_search = re.findall("href=\"\/drugs\/(DB\d+)#",webContentString)


    print(str(len(Drug_search)) + " Drugs")
    print(Drug_search)






    # For each of the Drugs check if a pdb-file exists
    # That means for each pdb open the web-page

    for i in range(0,len(Drug_search)):
        drug_url = 'https://www.drugbank.ca/drugs/' + Drug_search[i]

        path = "Drugs/"
        drug_url_file = path + Drug_search[i] + '.html'


        create_dir(path)

        if not os.path.isfile(drug_url_file):
            print("fetching data from " + drug_url + " ...")
            response = urllib2.urlopen(drug_url)
            webContent = response.read()
            f = open(drug_url_file, "w")
            f.write(webContent)
            f.close()
            

        webContent = open(drug_url_file, "r")
        webContentString = webContent.read()

        #print(webContentString)

        PDB_Entries_search = re.search(r"PDB Entries",webContentString)

        if PDB_Entries_search is not None:
            print("------------------------------------------------------------")
            print(Drug_search[i])
            #print(PDB_Entries_search)

            # extract the first pdb entry
            pdb_url_search = re.search(">PDB Entries.*?href=\"http:\/\/www\.rcsb\.org\/pdb\/explore\.do\?structureId=(.*?)\">.*?<\/a><\/span><span class='list-",webContentString)

            if pdb_url_search is not None:
                pdb_id = pdb_url_search.group(1)
                print(pdb_id)
                pdb_entry = "http://www.rcsb.org/pdb/explore.do?structureId=" + pdb_id
                print(pdb_entry)
                
                pdb_folder = target_name + "_pdbs/"
                create_dir(pdb_folder)

                pdb_download = "https://files.rcsb.org/view/" + pdb_id.capitalize() + ".pdb"
                pdb_file_name = pdb_folder + "/" + pdb_id.capitalize() + ".pdb"

                if not os.path.isfile(pdb_file_name):
                    response = urllib2.urlopen(pdb_download)
                    webContent = response.read()

                    f = open(pdb_file_name, "w")
                    f.write(webContent)
                    f.close()

#--------------------------------------------------------------------------------------------

target = 'P06401'
get_all_drugs_for_target(target)

target2 = 'P03372'
get_all_drugs_for_target(target2)

target3 = 'P10275'
get_all_drugs_for_target(target3)

target4 = 'P11511'
get_all_drugs_for_target(target4)


# Now merge all pdbs in one folder and create a file with the functionality
dirs= [name for name in os.listdir(".") if os.path.isdir(name)]

#print(dirs)
pdb_dirs = []
for i in range(0,len(dirs)):
    #print(dirs[i])
    if len(dirs[i].split("_pdbs")) == 2:
        #print("y")
        pdb_dirs.append(dirs[i])


print(pdb_dirs)

create_dir("merged")

for i in range(0,len(pdb_dirs)):
    print(pdb_dirs[i])
    fromDirectory = "./" + pdb_dirs[i] + "/"
    toDirectory = "./merged/"

    copy_tree(fromDirectory, toDirectory)


# create a file with the functionality
merged_files = get_files_in_dir("./merged/")
print(len(merged_files))

out = dict({}) 
for i in range(0,len(pdb_dirs)):
    pdb_files = get_files_in_dir("./" + pdb_dirs[i])
    print(pdb_files)
    for j in range(0,len(pdb_files)):
        if pdb_files[j] not in out:
            out[os.path.splitext(pdb_files[j])[0]] = pdb_dirs[i]
            print(os.path.splitext(pdb_dirs[i])[0], "_")
        else:
            print("functionality intersects! Skipping this protein!")


print(out)
print(len(merged_files))
print(len(out))




fout = "labels.txt"
fo = open(fout, "w")
fo.write("\"name\" \"label\"\n")
for k, v in out.items():
    fo.write("\"" + str(k) + "\"" + ' ' + "\"" + str(v) + "\"" + '\n')

fo.close()








