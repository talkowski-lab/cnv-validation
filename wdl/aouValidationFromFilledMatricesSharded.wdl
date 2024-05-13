version 1.0

import "genomeStripIRS.wdl" as gsirs
import "Structs.wdl"

workflow aouArrayValidation {

    input {
        Array[String] samples
        Array[File] filled_lrr_by_contig
        Array[File] per_contig_subset_gatk_sv_vcf
        Array[File] per_contig_subset_gatk_sv_vcf_idx
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

    Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]

    scatter (i in range(length(contigs))) {

        call ScatterVcf {
            input:
                vcf=per_contig_subset_gatk_sv_vcf[i],
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
                    array=filled_lrr_by_contig[i],
                    gs_tarball=gs_tarball,
                    array_validation_docker=array_validation_docker,
                    runtime_attr_override = runtime_attr_genome_strip_irs
            }
        }
    }

    call concatIrsReports {
        input:
            reports=flatten(genomeStripIRS.report),
            prefix=prefix,
            array_validation_docker=array_validation_docker,
            runtime_attr_override=runtime_attr_concat_irs_reports
    }

    output {
        Array[File] irs_vcf = flatten(genomeStripIRS.vcf)
        File irs_report = concatIrsReports.concat_report
    }

}


task concatIrsReports {
    input {
        Array[File] reports
        String prefix
        String array_validation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 64,
        disk_gb: 80,
        boot_disk_gb: 30,
        preemptible: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File concat_report = "~{prefix}.all_contigs.irs.report.tsv.gz"
    }

    command <<<
        set -eu

        echo "Merging IRS Report files"

        zcat ~{reports[0]} | head -n1 > ~{prefix}.all_contigs.irs.report.tsv

        set -o pipefail  # set after head -n1

        # assume already sorted by and within contig
        while read SHARD; do
            zcat $SHARD | tail -n+2 >> ~{prefix}.all_contigs.irs.report.tsv
        done < ~{write_lines(reports)}
        bgzip ~{prefix}.all_contigs.irs.report.tsv
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: array_validation_docker
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
