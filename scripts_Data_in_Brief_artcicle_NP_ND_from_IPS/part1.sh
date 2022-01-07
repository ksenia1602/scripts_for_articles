fastq="/data5/bio/runs-anufrieva/data/MA/IMG/178.213.246.140/root/projects/IMG/geo_neuron1"
trim="/data5/bio/runs-anufrieva/data/MA/IMG/trimming"
trimmomatic_prog="/data5/bio/runs-anufrieva/tools/Trimmomatic-0.35/trimmomatic-0.35.jar"
annotation="/data5/bio/runs-anufrieva/genome/salmon_hg38/salmon_index2"
quants="/data5/bio/runs-anufrieva/data/MA/IMG/salmon2"
salmon_prog="/data5/bio/runs-anufrieva/tools/salmon-0.12.0_linux_x86_64/bin/"

for j in {1..18}
do
#Trimming
      cd ${fastq}
      file1=($(find  ${fastq}/${j}_S${j}_R1_001.fastq.gz ))
      file2=($(find  ${fastq}/${j}_S${j}_R2_001.fastq.gz ))
      java -jar ${trimmomatic_prog} \
      PE -threads 20 -phred33 ${file1} ${file2} \
      ${trim}/IMG${j}_R1.trim.fastq.gz  R1_union_unpaired3.fastq.gz \
      ${trim}/IMG${j}_R2.trim.fastq.gz R2_union_unpaired3.fastq.gz \
      HEADCROP:1 \
      CROP:74 
      
#quasi-mapping
      ${salmon_prog}/salmon quant  \
      --libType A \
      -1 ${trim}/IMG${j}_R1.trim.fastq.gz  \
      -2 ${trim}/IMG${j}_R2.trim.fastq.gz  \
      --validateMappings \
      -p 20 \
      --index ${annotation} \
      -o ${quants}/${j}
done
