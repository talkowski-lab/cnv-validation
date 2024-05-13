version 1.0

import "genomeStripIRS.wdl" as gsirs
import "Structs.wdl"

workflow GenomeStripIRSSharded {

    input {
        Array[String] samples
        File filled_lrr
        File per_contig_subset_gatk_sv_vcf
        String prefix

        File primary_contigs_fai
        File genome
        File genome_index
        File genome_dict

        String scripts
        File gs_tarball
        String array_validation_docker

        Int records_per_shard
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_override_scatter
        RuntimeAttr? runtime_attr_genome_strip_irs
        RuntimeAttr? runtime_attr_concat_irs_reports
    }

    call ScatterVcf {
        input:
            vcf=per_contig_subset_gatk_sv_vcf,
            records_per_shard = records_per_shard,
            prefix = "~{prefix}.scatter_vcf",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_scatter
    }

    scatter (j in range(length(ScatterVcf.shards))) {
        call gsirs.genomeStripIRS {
            input:
                input_file=ScatterVcf.shards[j],
                prefix=prefix,
                genome=genome,
                genome_index=genome_index,
                genome_dict=genome_dict,
                array=filled_lrr,
                gs_tarball=gs_tarball,
                array_validation_docker=array_validation_docker,
                runtime_attr_override = runtime_attr_genome_strip_irs
        }
    }

    output {
        Array[File] vcf = genomeStripIRS.vcf
        Array[File] report = genomeStripIRS.report
    }

}


task ScatterVcf {
    input {
        File vcf
        String prefix
        Int records_per_shard
        Int? threads = 1
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 2,
                                      preemptible: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu, runtime_default.cpu])
        preemptible: select_first([runtime_override.preemptible, runtime_default.preemptible])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        # in case the file is empty create an empty shard
        bcftools view -h ~{vcf} | bgzip -c > ~{prefix}.0.vcf.gz
        bcftools +scatter ~{vcf} -o . -O z -p ~{prefix}. --threads ~{threads} -n ~{records_per_shard}

        ls ~{prefix}.*.vcf.gz | sort -k1,1V > vcfs.list
        i=0
        while read vcf; do
            shard_no=`printf %06d $i`
            mv ${vcf} ~{prefix}.shard_${shard_no}.vcf.gz
            tabix -p vcf ~{prefix}.shard_${shard_no}.vcf.gz
            i=$((i+1))
        done < vcfs.list
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.gz")
        Array[File] shards_idx = glob("~{prefix}.shard_*.vcf.gz.tbi")
    }
}
