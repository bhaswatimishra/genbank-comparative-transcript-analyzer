from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
import pandas as pd
import subprocess
from Bio import Phylo


def genbank_accession(gene,correct_organism) :
    try :
          Entrez.email = "bhaswatimishra154@gmail.com"

          handle = Entrez.esearch(
                   db="nucleotide",
                   term=f"{gene}[Gene] AND {correct_organism}[Organism] AND MANE Select[keyword]"
                   )

          record = Entrez.read(handle)
          if not record['IdList'] :
                print("No ids found for the desired gene.")
                return False
          elif len(record['IdList']) >1 :
                print("more than 1 ids found")
                print(record['IdList'])
                return False
          else :
               url = Entrez.esummary(
                     db = "nucleotide",
                     id = f"{record['IdList'][0]}",
                     rettype = "xml"
                     )
               new_record = Entrez.read(url)
               nt_id = new_record[0]['AccessionVersion']
               nt_url = Entrez.efetch(
                    db = "nucleotide",
                    id = nt_id,
                    rettype = "gb",
                    retmode = "text"
                    )
               nt_record = nt_url.read()
               with open(f"{gene}_{correct_organism}.gb", "w") as f :
                       f.write(nt_record)
               return True
    except Exception as E :
         print (f"Error : {E}")

def genbank_datafetch(gene,correct_organism,protein_id_list) :
     try :
          if os.path.exists(f"{gene}_{correct_organism}.gb") :
               comparision[correct_organism] = {}
               comparision[correct_organism]['protein_data'] = {}
               record = SeqIO.read(f"{gene}_{correct_organism}.gb","genbank")
               comparision[correct_organism]['seq_len'] = len(record.seq)
               comparision[correct_organism]['GC_proportion'] = round(gc_fraction(record.seq) * 100,2)
               comparision[correct_organism]['topology'] = record.annotations.get('topology')
               count = 0
               for feature in record.features :
                    if feature.type == "exon" :
                         count += 1
                    elif feature.type == "CDS" :
                         protein_id = feature.qualifiers['protein_id'][0]
                         comparision[correct_organism]['protein_data'] = {}
                         comparision[correct_organism]['protein_data']['protein_id'] = protein_id
                         translation = feature.qualifiers['translation'][0]
                         protein = ProteinAnalysis(translation)
                         comparision[correct_organism]['protein_data']['mol_weight'] = round(protein.molecular_weight(),2)
                         comparision[correct_organism]['protein_data']['total_AA'] = protein.length
                         comparision[correct_organism]['protein_data']['isoelectric_point'] = round(protein.isoelectric_point(),2)
                         comparision[correct_organism]['protein_data']['aromaticity'] = round(protein.aromaticity(),2)
                         protein_id_list.append(protein_id)
                    elif feature.type == "source" :
                         comparision[correct_organism]['chromosome'] = feature.qualifiers.get("chromosome", ["Unknown"])[0]
                    else :
                         continue
               if count == 0 :
                   comparision[correct_organism]['exon_count'] = "N/A"
               else :
                    comparision[correct_organism]['exon_count'] = count 
               return True
     except Exception as E :
          print(f" Error 2 : {E}")

def genbank_to_fasta(gene,correct_organism) :
     if os.path.exists(f"{gene}_{correct_organism}.gb") :
          try :
               SeqIO.convert( f"{gene}_{correct_organism}.gb", "genbank" , f"{gene}_{correct_organism}.fasta", "fasta")
               print( f"Fasta file created successfully")
               print(f"file name : {gene}_{correct_organism}.fasta")
          except Exception as E :
               print("Couldn't created fasta file")
               print(f"Error : {E}")

def blast_comparision(gene,correct_organism,new_id) :
  try :
     nt_url = Entrez.efetch(
                    db = "nucleotide",
                    id = new_id,
                    rettype = "gb",
                    retmode = "text"
                    )
     nt_record = nt_url.read()
     with open(f"{gene}_{correct_organism}.gb", "w") as f :
          f.write(nt_record)
     try :
         genbank_datafetch(gene,correct_organism,protein_id_list)
     except Exception as E :
          print(f"Error : {E}")
          print("failed to load genbank analysis")
  except Exception as E :
       print("invalid id")
       print(f"error : {E}")


def comparision_report(comparision) :
     report = {}
     report['features'] = [ 'length',
                           'GC count',
                           'topology',
                           'protein id',
                           'mol weight of protein',
                           'total AA',
                           'isoelectric point',
                           'aromaticity',
                           'chromosome no',
                           'exon count']
     comp_list = list(comparision.keys())
     for i in range(0,len(comp_list)) :
         index = comp_list[i] 
         report[index] = []
         report[index].append(comparision[index]['seq_len'])
         report[index].append(comparision[index]['GC_proportion']) 
         report[index].append(comparision[index]['topology'])
         report[index].append(comparision[index]['protein_data']['protein_id']) 
         report[index].append(comparision[index]['protein_data']['mol_weight']) 
         report[index].append(comparision[index]['protein_data']['total_AA']) 
         report[index].append(comparision[index]['protein_data']['isoelectric_point']) 
         report[index].append(comparision[index]['protein_data']['aromaticity'])
         report[index].append(comparision[index]['chromosome']) 
         report[index].append(comparision[index]['exon_count']) 
     df = pd.DataFrame(report)
     df.to_csv("transcript_comparator.csv", index=False)
     print(f"Report saved successfully as : transcript_comparator.csv")

def genbank_to_allignfasta(protein_id_list) :
     for protein_id in protein_id_list :
          try :
             handle = Entrez.efetch(
                   db = "protein",
                   id = f"{protein_id}",
                   rettype = "fasta",
                   retmode = "text"
                   )
             record = handle.read()
             with open("all_sample.fasta", "a") as f :
               f.write(record)
          except Exception as E :
               print("couldn't convert fasta file for allignment")
               print(E)

def make_alignment(fasta_file) :
     try :
        if os.path.exists(fasta_file) :
          subprocess.run(
               ['muscle',
                '-align',
                fasta_file,
                '-output',
                'aligned.fasta']
          )
     except Exception as E :
          print("failed to create alignment file")
          print(E)

def phylogenic_tree(aligner_fasta) :
    try:
        if os.path.exists(aligner_fasta):
            with open("tree.nwk", "w") as tree:
                subprocess.run(
                    ["fasttree", aligner_fasta],
                    stdout = tree ,
                    check=True
                )

            print("Phylogenetic tree created successfully.")
            print("File name : tree.nwk")

        else:
            print("Alignment file not found.")

    except Exception as E:
        print("Failed to build phylogenetic tree.")
        print(E)
        

def visualize_tree(tree_file):
    try:
        tree = Phylo.read(tree_file, "newick")
    
        Phylo.draw_ascii(tree, column_width= 80 )
       
    except Exception as E:
        print("Unable to display tree")
        print(E)

comparision = {}
protein_id_list = []
if os.path.exists("all_sample.fasta"):
    os.remove("all_sample.fasta")
gene = input("Enter the gene name :").upper()
correct_organism1 = "Homo sapiens"
print("By default the organism name is set as human")
if genbank_accession(gene, correct_organism1) :
    genbank_datafetch(gene, correct_organism1,protein_id_list)
    genbank_to_fasta(gene, correct_organism1)
    
else:
    exit()
print("use the fasta file to run BLAST and from the matching sequences enter the desired accession id")
new_id = input("Enter the id :").upper()
organism_name = input("Enter the name of the organism :").strip()
correct_name = organism_name.replace(" ","_")
blast_comparision(gene,correct_name,new_id)

while True :
     more = input("more species id to add ?").upper()
     if more == "YES" :
          new_id = input("Enter the id :").upper()
          organism_name = input("Enter the name of the organism :").strip()
          correct_name = organism_name.replace(" ","_")
          blast_comparision(gene,correct_name,new_id)
     else :
          comparision_report(comparision)
          genbank_to_allignfasta(protein_id_list)
          make_alignment("all_sample.fasta")
          phylogenic_tree("aligned.fasta")
          visualize_tree("tree.nwk")
          break
