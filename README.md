# ncbi_metadata_collection
Collection and organization of metadata associated with NCBI nucleotide and genome sequences for bioinformatics and comparative genomic analyses, including BEAST, Phylogenetic analysis, and geographic distribution studies
#Install the required dependencies:

conda install -c conda-forge ncbi-datasets-cli biopython pandas requests -y

#Clone the GitHub repository:

git clone https://github.com/ktronimia/ncbi-metadata-collection.git

#Navigate to the repository directory:

cd ncbi-metadata-collection

#Run the Python script:

python ncbi_metadata_collection_01.py

python ncbi_metadata_collection_universal.py

python metadata_collection_using_accession.py (note: before run "metadata_collection_using_accession.py" , create accession_list.txt, format mentioned inside the code)
