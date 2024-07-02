nextflow.enable.dsl=2
//include { bowtie2_align_PE_S_wf } from '/home/tfallon/source/nextflow/bowtie2/bowtie2.nf'

process cat_collected_files {
input:
 path(files)
output:
 path("cat-ed_*")
shell:
'''
##filename=$(basename -- !{files[0]})
##extension="${filename##*.}"
##filename="${filename%.*}"

cat !{files} > cat-ed_!{files[0]}
'''
}

process seqkit_sum_rename_fasta {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
scratch 'ram-disk'
conda 'seqkit'
input:
 path(fasta)
output:
 path("seqkit.*")  
tag "$fasta"
shell:
'''
SUM=$(seqkit sum !{fasta} | cut -f 1)
seqkit replace -p "^" -r "${SUM}__" !{fasta} > ${SUM}__!{fasta}
'''
}

process trf_checksum_fasta {
executor 'local'
conda 'seqkit openssl coreutils'
publishDir "results",pattern:"input*.checksum.txt",mode:"copy",overwrite:"true"
input:
 path genome
output:
 tuple env(FASCHK),path("*-*-*-*__${genome}"), emit:checksum_and_genome
 path "input.*.checksum.txt"
tag "${genome}"
shell:
'''
### Recommend use seqkit sum instead
f=!{genome}
FCHK=$(cat $f | openssl md5 | cut -f 2 -d " " | cut -c1-4) ##First 4 characters of a file contents md5 checksum
IDCHK=$(seqkit seq -n -i $f | sort | openssl md5 | cut -f 2 -d " " | cut -c1-6) ##First 6 characters of a md5 checksum of the sorted, concatenated FASTA IDs
SEQCHK=$(seqkit seq -u $f | seqkit sort -s | seqkit seq -s | openssl md5 | cut -f 2 -d " " | cut -c1-6) ##First 6 characters of a md5 checksum of the sorted, concatenated, uppercase, FASTA sequence
ESEQCHK=$(seqkit sort -s $f | seqkit seq -s | openssl md5 | cut -f 2 -d " " | cut -c1-4) ##First 4 characters of a md5 checksum of the sorted, concatenated FASTA sequence
FASCHK="${FCHK}-${IDCHK}-${SEQCHK}-${ESEQCHK}"
ln -s $f ${FASCHK}__${f}
echo "faschk:${FASCHK}__${f}"
echo "faschk:${FASCHK}__${f}" > input.${f}.checksum.txt
'''
}

process gt_tidysort {
input:
 path gff_file
output:
 path "results/${gff_file}"
shell:
'''
mkdir results
gt gff3 -tidy -sort -retainids !{gff_file} > results/!{gff_file}
'''
}

process detab {
input:
 path toDetab
output:
 path "done/${toDetab}"
shell:
"""
mkdir done
sed "s^\t^ ^g" !{toDetab} > done/!{toDetab}
"""
}


process debar {
input:
 path toDebar
output:
 path "done/${toDebar}"
shell:
"""
mkdir done
sed 's^|^_^g' !{toDebar} > done/!{toDebar}
"""
}

process decolon {
input:
 path toDecolon
output:
 path "done/${toDecolon}"
shell:
"""
mkdir done
sed 's^:^=^g' !{toDecolon} > done/!{toDecolon}
"""
}

process dummy_publish_path {
publishDir "results/", mode: 'link',overwrite:'true'
input:
 path theDummy
output
 path "${theDummy}"
shell:
"""
echo "mainly just here to publish the file"
"""
}

workflow rename_and_map_PE_S_wf {

fasta = Channel.fromPath(params.fasta)
forward = Channel.fromPath(params.forward)
reverse = Channel.fromPath(params.reverse)
singles = Channel.fromPath(params.singles)

pairedReads = forward.combine(reverse)
singles_merged = singles.collectFile(name:'merged_singles.fastq.gz')

seqkit_sum_rename_fasta(fasta)

//bowtie2_align_PE_S_wf(pairedReads,singles_merged,seqkit_sum_rename_fasta.out)

}
