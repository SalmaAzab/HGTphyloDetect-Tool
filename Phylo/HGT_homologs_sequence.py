import re
import sys
import os
from Bio import SeqIO
from Bio import Entrez
from ete3 import NCBITaxa

def parse_NCBI(filename):
    with open(filename, "r") as infile:
        lines = infile.readlines()

    accession_number = list()
    accession_similarity = dict()
    for line in lines:
        accession = line.strip("\n").split("\t")[2]
        similarity = line.strip("\n").split("\t")[3]
        accession_number.append(accession)
        accession_similarity[accession] = float(similarity)

    return accession_number, accession_similarity

def get_refSeq(gene):
    # get the related protein sequence according to protein identifier
    with open(f"./input/{gene}.fasta", "r") as handleGene:
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta"):
            proteinSeq[record.id] = str(record.seq)

    return proteinSeq

def getTaxid(accession):
    Entrez.email = "abcd@ncbi.org"
    handle = Entrez.efetch(db='protein', id=accession, rettype='gb')
    record = SeqIO.read(handle, 'genbank')
    if record.features[0].qualifiers['db_xref'][0].split(":")[0] == 'taxon':
        taxid = record.features[0].qualifiers['db_xref'][0].split(":")[1]
        organism = record.features[0].qualifiers['organism'][0]

    return taxid

def getSeq(accession):
    Entrez.email = "abcd@ncbi.org"
    handle = Entrez.efetch(db='protein', id=accession, rettype='gb')
    record = SeqIO.read(handle, 'genbank')
    seq = record.seq

    return seq

def main():
    filename = sys.argv[1]
    gene = os.path.splitext(os.path.basename(filename))[0]

    print('This is gene', gene + ' ---------------')
    if os.path.exists(f"./input/{gene}.txt"):
        print(f"File ./input/{gene}.txt exists.")
        accession_number, accession_similarity = parse_NCBI(f"./input/{gene}.txt")
        print(f"Parsed {len(accession_number)} accession numbers.")
    else:
        print(f"Warning: please run BLASTP first! File ./input/{gene}.txt not found.")
        return  # Exit the function if the file does not exist

    ncbi = NCBITaxa()
    id_seq = dict()
    for accession in accession_number:
        try:
            seq = str(getSeq(accession))
        except:
            continue

        taxid = getTaxid(accession)
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())

        try:
            taxid2name = ncbi.get_taxid_translator([ranks2lineage['phylum'], ranks2lineage['species']])
            if taxid2name[ranks2lineage["phylum"]] != "Ascomycota":
                accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@" + taxid2name[ranks2lineage["phylum"]]
                print(accession)
                id_seq[accession] = seq
                continue
        except:
            pass

        try:
            taxid2name = ncbi.get_taxid_translator([ranks2lineage['phylum'], ranks2lineage['subphylum'], ranks2lineage['species']])
            if accession_similarity[accession] < 80:
                if taxid2name[ranks2lineage["subphylum"]] == "Saccharomycotina":
                    accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@Saccharomycotina"
                if taxid2name[ranks2lineage["subphylum"]] != "Saccharomycotina" and taxid2name[ranks2lineage["phylum"]] == "Ascomycota":
                    accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@other_Ascomycota"
                if taxid2name[ranks2lineage["phylum"]] != "Ascomycota":
                    accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@" + taxid2name[ranks2lineage["phylum"]]
                print(accession)
                id_seq[accession] = seq
            else:
                continue
        except:
            continue

    gene_seq = get_refSeq(gene)
    gene_query = "QUERY_" + 'Saccharomyces_cerevisiae_' + gene
    print(gene_query)
    print(gene_seq)
    id_seq[gene_query] = gene_seq

    with open(f"./input/{gene}_homologs.fasta", "w") as outfile:
        for accession, seq in id_seq.items():
            outfile.write(f">{accession}\n")
            outfile.write(f"{seq}\n")

if __name__ == "__main__":
    main()
