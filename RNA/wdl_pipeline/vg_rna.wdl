#version draft-2
# Apparently cromwell 36.1 doesn't parse the version line in wdl files if the wdl version is a draft version.
# Cromwell 36.1 wdl parser defaults to draft-2 if a `version` line is not given.


workflow vgRNAPipeline {
    File INPUT_FASTQ_READ_PAIR1
    File INPUT_FASTQ_READ_PAIR2
    File INPUT_REF_FASTA
    File INPUT_REF_FASTA_INDEX
    File INPUT_SAMPLE_PHASED_VCF
    File INPUT_SAMPLE_PHASED_VCF_INDEX
    File INPUT_GENE_GTF
    String CONTIG_NAME
    String SAMPLE_NAME
    String VG_CONTAINER
    
    call constructVGandGBWT {
        input:
            in_ref_fasta=INPUT_REF_FASTA,
            in_ref_fasta_index=INPUT_REF_FASTA_INDEX,
            in_sample_phased_vcf=INPUT_SAMPLE_PHASED_VCF,
            in_sample_phased_vcf_index=INPUT_SAMPLE_PHASED_VCF_INDEX,
            in_contig_name=CONTIG_NAME,
            in_sample_name=SAMPLE_NAME,
            in_vg_container=VG_CONTAINER
    }
    call constructGeneGraphs {
        input:
            in_vg_file=constructVGandGBWT.output_vg_file,
            in_gbwt_file=constructVGandGBWT.output_gbwt_file,
            in_gene_gtf=INPUT_GENE_GTF,
            in_contig_name=CONTIG_NAME,
            in_sample_name=SAMPLE_NAME,
            in_vg_container=VG_CONTAINER
    } 
    call runVGMPMAP {
        input:
            in_fastq_read_pair1=INPUT_FASTQ_READ_PAIR1,
            in_fastq_read_pair2=INPUT_FASTQ_READ_PAIR2,
            in_xg_file=constructGeneGraphs.output_xg_file,
            in_gcsa_files=constructGeneGraphs.output_gcsa_files,
            in_dist_file=constructGeneGraphs.output_dist_file,
            in_snarls_file=constructGeneGraphs.output_snarls_file,
            in_sample_name=SAMPLE_NAME,
            in_vg_container=VG_CONTAINER
    } 
}

task constructVGandGBWT {
    File in_ref_fasta
    File in_ref_fasta_index
    File in_sample_phased_vcf
    File in_sample_phased_vcf_index
    String in_contig_name
    String in_sample_name
    String in_vg_container

    String dollar = "$" 
    command <<< 
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        vg construct -a -r ${in_ref_fasta} -v ${in_sample_phased_vcf} -R ${in_contig_name} -t 16 -p > ${in_sample_name}.vg && \
        mkdir vg_index_tmp && \
        vg index -t 16 -p -b vg_index_tmp -v ${in_sample_phased_vcf} --gbwt-name ${in_sample_name}.gbwt ${in_sample_name}.vg
    >>> 
    output {
        File output_vg_file = glob("${in_sample_name}.vg")[0]
        File output_gbwt_file = glob("${in_sample_name}.gbwt")[0]
    }
    runtime {
        docker: in_vg_container
    }
}

task constructGeneGraphs {
    File in_vg_file
    File in_gbwt_file
    File in_gene_gtf
    String in_contig_name
    String in_sample_name
    String in_vg_container

    String dollar = "$" 
    command <<< 
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        vg mod -t 16 -r ${in_contig_name} ${in_vg_file} > ${in_sample_name}_tmp.vg && \
        mv ${in_sample_name}_tmp.vg ${in_sample_name}.vg && \
        sed -i 's/chr//g' ${in_gene_gtf} && \
        vg rna -p -t 16 -c -a -n ${in_gene_gtf} -l ${in_gbwt_file} -b ${in_sample_name}_gencode.gbwt -g ${in_sample_name}_gencode.gam -f ${in_sample_name}_gencode.fa ${in_sample_name}.vg | vg mod -t 16 -r ${in_contig_name} -I - | vg mod -t 16 -N - | vg ids -c -s - > ${in_sample_name}_gencode_paths_genes_sort.vg && \
        vg snarls -t ${in_sample_name}_gencode_paths_genes_sort.vg > ${in_sample_name}_gencode_paths_genes_sort.snarls && \
        mkdir vg_index_tmp && \
        vg index -t 16 -p -b vg_index_tmp --xg-name ${in_sample_name}_gencode_paths_genes_sort.xg --gcsa-out ${in_sample_name}_gencode_paths_genes_sort.gcsa -s ${in_sample_name}_gencode_paths_genes_sort.snarls -w 100 -j ${in_sample_name}_gencode_paths_genes_sort.dist ${in_sample_name}_gencode_paths_genes_sort.vg
    >>> 
    output {
        File output_xg_file = glob("${in_sample_name}_gencode_paths_genes_sort.xg")[0]
        Array[File] output_gcsa_files = glob("${in_sample_name}_gencode_paths_genes_sort.gcsa*")
        File output_dist_file = glob("${in_sample_name}_gencode_paths_genes_sort.dist")[0]
        File output_snarls_file = glob("${in_sample_name}_gencode_paths_genes_sort.snarls")[0]
    }
    runtime {
        docker: in_vg_container
    }
}

task runVGMPMAP {
    File in_fastq_read_pair1
    File in_fastq_read_pair2
    File in_xg_file
    Array[File] in_gcsa_files
    File in_dist_file
    File in_snarls_file
    String in_sample_name
    String in_vg_container
    
    File in_gcsa_file = in_gcsa_files[0]
    File in_gcsa_lcp_file = in_gcsa_files[1]
    String dollar = "$" 
    command <<< 
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
       
        time vg mpmap \
            -t 16 \
            --fastq ${in_fastq_read_pair1} \
            --fastq ${in_fastq_read_pair2} \
            --xg-name ${in_xg_file} \
            --gcsa-name ${in_gcsa_file} \
            --dist-name ${in_dist_file} \
            --snarls ${in_snarls_file} \
            --min-dist-cluster \
            > ${in_sample_name}_gencode_paths_genes_sort.gamp && \
        vg view -K -G ${in_sample_name}_gencode_paths_genes_sort.gamp > ${in_sample_name}_gencode_paths_genes_sort.gam
    >>> 
    output {
        File output_gamp_file = glob("${in_sample_name}_gencode_paths_genes_sort.gamp")[0]
        File output_gam_file = glob("${in_sample_name}_gencode_paths_genes_sort.gam")[0]
    }
    runtime {
        docker: in_vg_container
    }
}

