version 1.0

workflow pipeline_1 {

    input {
        Array[String] Samples
        String file_path
        File ref_dir
        String fasta
        String SENTIEON_LICENSE
        String group
        String platform
        String? duplex_umi
        String? read_structure
        Boolean qc=false
        File? fastq_1
        File? fastq_2
        
    }

    scatter (sample_id in Samples){
        if (qc){
              call fastp{
                input:
                sample_id = sample_id,
                file_path = file_path,
                docker = "registry.cn-shanghai.aliyuncs.com/hcc1395_aliyun/fastp:0.19.6"
            }
              call mapping as mapping1{
                input:
                ref_dir = ref_dir,
                fasta = fasta,
                fastq_1 = fastp.out1,
                fastq_2 = fastp.out2,
                SENTIEON_LICENSE = SENTIEON_LICENSE,
                group = group,
                sample = sample_id,
                platform = platform,
                duplex_umi = duplex_umi,
                read_structure = read_structure,
                docker = "registry.cn-shanghai.aliyuncs.com/hcc1395_aliyun/sentieon-genomics:v202112.05",

            }
        }

        if (fastq_1 != ""){
              call mapping as mapping2{
                input:
                ref_dir = ref_dir,
                fasta = fasta,
                fastq_1 = fastq_1,
                fastq_2 = fastq_2,
                SENTIEON_LICENSE = SENTIEON_LICENSE,
                group = group,
                sample = sample_id,
                platform = platform,
                duplex_umi = duplex_umi,
                read_structure = read_structure,
                docker = "registry.cn-shanghai.aliyuncs.com/hcc1395_aliyun/sentieon-genomics:v202112.05",

            }

        }


    }



}

task fastp {

    input {
        
        # I/O options
        String sample_id
        String file_path

        Boolean? phred64 = false 
        Boolean? fix_mgi_id = false

        String? adapter_sequence
        String? adapter_sequence_r2

        Int? reads_to_process # specify how many reads/pairs to be processed. Default 0 means process all reads.


        # excute env
        Int cpu = 2
        String docker
        String memory = "4G"
        String disks = "local-disk 50 cloud_ssd"

    }

    File in1 = file_path + sample_id +'_reads_1.fastq'
    File in2 = file_path + sample_id +'_reads_2.fastq'
    
    String out1_name = sample_id + "-clean" +'_reads_1.fastq' 
    String out2_name = sample_id + "-clean" +'_reads_2.fastq' 
    
    # reporting options
    String json = sample_id + "_fastp.json"
    String html = sample_id + "_fastp.html"
    String report_title = sample_id + " \'fastp report\'"

    command <<<

        # basic command
        /opt/conda/bin/fastp \
        --in1 ~{in1} \
        --in2 ~{in2} \
        --out1 ~{out1_name} \
        --out2 ~{out2_name} \
        --json ~{json} \
        --html ~{html} \
        --report_title ~{report_title} \
        
        # options 
        ~{ true="--phred64 " false="" phred64 } \
        ~{ "--reads_to_process " + reads_to_process } \
        ~{ true="--fix_mgi_id " false="" fix_mgi_id } \
        ~{ "--adapter_sequence " + adapter_sequence } \
        ~{ "--adapter_sequence_r2 " + adapter_sequence_r2 }

    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker
    }

    output {
        File out1 = out1_name
        File out2 = out2_name
        File json_report = json
        File html_report = html
    }

}


task mapping {
    input {
          File ref_dir
          String fasta
          File? fastq_1
          File? fastq_2

          String SENTIEON_LICENSE
          String group
          String sample
          String platform
          String? duplex_umi
          String? read_structure
          Int cpu = 2
          String docker
          String memory = "4G"
          String disks = "local-disk 50 cloud_ssd"

    }


  command <<<
    set -o pipefail
    set -e
    export SENTIEON_LICENSE=${SENTIEON_LICENSE}
    nt=$(nproc)
    
    if [ ~{read_structure} ]; then
      if [ ~{duplex_umi} == "true" ]; then
        READ_STRUCTURE="-d ~{read_structure}"
      fi
      sentieon umi extract $READ_STRUCTURE ~{fastq_1} ~{fastq_2} | \
      sentieon bwa mem -p -C -R "@RG\tID:~{group}\tSM:~{sample}\tPL:~{platform}" -t $nt -K 10000000 ~{ref_dir}/~{fasta} - | \
      sentieon umi consensus -o ~{sample}.umi_consensus.fastq.gz
      
      sentieon bwa mem -p -C -R "@RG\tID:~{group}\tSM:~{sample}\tPL:~{platform}" -t $nt -K 10000000 ~{fasta} ~{sample}.umi_consensus.fastq.gz | \
      sentieon util sort --umi_post_process --sam2bam -i - -o ~{sample}.sorted.bam
    else
      sentieon bwa mem -R "@RG\tID:~{group}\tSM:~{sample}\tPL:~{platform}" \
      -t $nt -K 10000000 ~{ref_dir}/~{fasta} ~{fastq_1} ~{fastq_2} | \
      sentieon util sort -o ~{sample}.sorted.bam -t $nt --sam2bam -i -
    fi
    
  >>>

  runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker
  }
  output {
    File sorted_bam = "~{sample}.sorted.bam"
    File sorted_bam_index = "~{sample}.sorted.bam.bai"
  }
}

task pcard_tools {
  input{}
  command<<<
  >>>
  output{}
}


