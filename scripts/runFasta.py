#!/usr/bin/env python
# coding: utf-8

# ### The purpose of this notebook is to composite several bioinformatics techniques and tools into one place. 
# 
# Currently -- 1) user inputs gene of interest and variants of interest
#              2) gathers database and generates table from chembl 
#              3) user inputs table id selection with correct gene of interest
#              4) grabs fasta file from uniprot
#              5) generates file named with uniprot ID_mutation (EX: P49768_T380M.fasta)
#                  IF: wild type at that position matches input wild type then a swap occurs
#                  ELSE: mutation swap does not occur and file is just a copy of wild type
#              
# Goal      -- 1) ~~Get request from Polyphen 2 site, open results, and print results. ~~
#              2) Generate graphs/tables based on analysis from (...)
#        
# Author: Luke Terwilliger 
# Version: 1.3
# Last update: 05/08/22
# ###
from IPython import get_ipython
import os
import pandas as pd
from chembl_webresource_client.new_client import new_client
import json
import requests
import os
import shutil

gene = input("Enter the gene of interest EX: PSEN2: ")

mutationSwap = input("Enter mutation swaps (original amino acid Position New amino Acid) \n[EX T380M,T380W,T350M ]: ")
mutationSwap = mutationSwap.upper()

# get_ipython().system('pip install chembl_webresource_client')


# Windows / Linux(not tested)
# !pip install pandas
# MAC OS ARM 
#get_ipython().system('pip3 install cython')
#get_ipython().system('OPENBLAS="$(brew --prefix openblas)" MACOSX_DEPLOYMENT_TARGET=11.1 pip3 install numpy --no-use-pep517')
#get_ipython().system('OPENBLAS="$(brew --prefix openblas)" MACOSX_DEPLOYMENT_TARGET=11.1 pip3 install pandas --no-use-pep517')


# Target search 
target = new_client.target
#target_query = target.search(gene)
target_query = target.filter(target_synonym__icontains=gene)
targets = pd.DataFrame.from_dict(target_query)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 3000)
pd.set_option('display.colheader_justify', 'center')
pd.set_option('display.precision', 3)
print(targets[['cross_references','organism','pref_name','target_chembl_id']])


# pick table id number from far left column and insert below
user_select = int(input("Please choose the table id(far left starting at 0) for your target gene: "))
selected_target = targets.target_chembl_id[user_select]
#print(selected_target)


# clean up cross references column to include only uniprot protein id
selected_target = targets.target_components[user_select]
selected_target[0]=json.dumps(selected_target)
selected_target=selected_target[0].split(":")
# chaos if uncommented, but possible other usable data here

#for i in range (0,len(selected_target)):
   #print(selected_target[i])
selected_target=selected_target[1].split(",")
selected_target = selected_target[0]
selected_target=selected_target.replace('"',"")
selected_target=selected_target.strip()
#print(selected_target)

webname = "https://www.uniprot.org/uniprot/" + selected_target
fasta = ".fasta"
fastaWeb = webname+fasta
#print(webname)
#print(fastaWeb)


# #Will open the webpage if needed, for testing
# import webbrowser
# 
# webbrowser.open(webname)

r = requests.get(fastaWeb, allow_redirects=True)

open(selected_target+fasta, 'wb').write(r.content)

mutationSwap=mutationSwap.replace(" ","")
#print(mutationSwap)
listMutations = mutationSwap.split(",")
listLen = len(listMutations)
#print(listMutations)


# save file as a new name 
# get source file name
src = selected_target+fasta
dest = ""
# get destination file name
for x in range (0,listLen):
    dest = selected_target+"_"+listMutations[x]+fasta
  #  print(src+" "+dest)
    shutil.copy(src, dest)

content=""
#print (listLen)
for x in range (0,listLen):
    #my_file = open(selected_target+"_"+listMutations[x]+fasta, "w")
   # print(selected_target+"_"+listMutations[x]+fasta)
    my_file = open(selected_target+"_"+listMutations[x]+fasta)
    content = my_file.read()
    my_file.close()
    

# gets the amino acid position from mutation list
for x in range (0,listLen):
    a_string = " ".join(listMutations[x])
    numbers = []
    for word in a_string.split():
        if word. isdigit():
            numbers.append(int(word))
    string_ints = [str(int) for int in numbers]
    str_of_ints = "". join(string_ints)
   # print(str_of_ints)

# this will loop through each mutation, write to a new file of the form (target_mutation.fasta), 
# swap the wild type amino acid with the mutant amino acid
for x in range(0,listLen):
    file1 = open(selected_target+fasta, 'r')
    file2 = open(selected_target+"_"+listMutations[x]+fasta, 'w')
    count = 0
    count2 = 0
    aminoAcids = "".join(listMutations[x])
    aminoWild = aminoAcids[0]
    aminoMutant = aminoAcids[len(aminoAcids)-1]
    a_string = " ".join(listMutations[x])
    numbers = []
    for word in a_string.split():
        if word. isdigit():
            numbers.append(int(word))
    string_ints = [str(int) for int in numbers]
    str_of_ints = "". join(string_ints)
    while True:
        count += 1
        # Get next line from file
        line = file1.readline()
        line=line.strip()
        listChar = list(line)
        # if line is empty
        # end of file is reached
        if not line:
            break
        if count > 1: 
            for i in range (0,len(listChar)):
                count2+=1
                
                if count2==int(str_of_ints):
                    if aminoWild == listChar[i]:
                        listChar[i]=aminoMutant 
                        print("Swapped HERE")
                    else:
                        print("NO SWAP, incorrect starting amino acid. File: " + selected_target+"_"+listMutations[x]+fasta)
                file2.write(listChar[i])
                if (i+1) % 60 == 0:
                    file2.write("\n")
        else:
            file2.write(line) 
            file2.write("\n")
    file1.close()
    file2.close()

# stuff that would be inserted into batch query form 
list = []
for x in range(0,listLen):
    aminoAcids = "".join(listMutations[x])
    aminoWild = aminoAcids[0]
    aminoMutant = aminoAcids[len(aminoAcids)-1]
    a_string = " ".join(listMutations[x])
    numbers = []
    for word in a_string.split():
        if word. isdigit():
            numbers.append(int(word))
    string_ints = [str(int) for int in numbers]
    str_of_ints = "". join(string_ints)
    list.append(selected_target + " "+str_of_ints +" "+ aminoWild + " "+aminoMutant)
    print(selected_target + " "+str_of_ints +" "+ aminoWild + " "+aminoMutant)

