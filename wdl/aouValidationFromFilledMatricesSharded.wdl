version 1.0

import "GenomeStripIRSSharded.wdl" as gsirs_sharded
import "Structs.wdl"

workflow aouArrayValidation {

    input {
        Array[String] samples
        Array[File] filled_lrr_by_contig
        Array[File] per_contig_gatk_sv_vcf
        Array[File] per_contig_gatk_sv_vcf_idx
        String prefix

        Int max_cnv_size=10000000
        String? min_cnv_size  # lower bound on SVLEN to evaluate against arrays. Default: 50000 (50kb)
        Int? max_ac  # maximum allele count to evaluate against arrays

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
        RuntimeAttr? runtime_attr_subset_gatk_sv
        RuntimeAttr? runtime_attr_genome_strip_irs
        RuntimeAttr? runtime_attr_concat_irs_reports
    }

    Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]

    scatter (i in range(length(contigs))) {

        call subsetGATKSV {
            input:
                gatk_sv_vcf=per_contig_gatk_sv_vcf[i],
                gatk_sv_vcf_idx=per_contig_gatk_sv_vcf_idx[i],
                sample_list=write_lines(samples),
                prefix=prefix,
                max_ac=max_ac,
                min_cnv_size=select_first([min_cnv_size, "50000"]),
                max_cnv_size=max_cnv_size,
                chromosome=contigs[i],
                array_validation_docker=array_validation_docker,
                scripts=scripts,
                runtime_attr_override = runtime_attr_subset_gatk_sv
        }

        call gsirs_sharded.GenomeStripIRSSharded as GenomeStripIRSSharded {
            input:
                per_contig_subset_gatk_sv_vcf=subsetGATKSV.subset_vcf,
                prefix=prefix,
                genome=genome,
                genome_index=genome_index,
                genome_dict=genome_dict,
                filled_lrr=filled_lrr_by_contig[i],
                gs_tarball=gs_tarball,
                records_per_shard=records_per_shard,
                scripts=scripts,
                array_validation_docker=array_validation_docker,
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override_scatter = runtime_attr_override_scatter,
                runtime_attr_genome_strip_irs = runtime_attr_genome_strip_irs
        }
    }

    call concatIrsReports {
        input:
            reports=flatten(GenomeStripIRSSharded.report),
            prefix=prefix,
            array_validation_docker=array_validation_docker,
            runtime_attr_override=runtime_attr_concat_irs_reports
    }

    output {
        Array[File] irs_vcf = flatten(GenomeStripIRSSharded.vcf)
        File irs_report = concatIrsReports.concat_report
    }

}

task RemoveVeryLargeEvents {
    input {
        File gatk_sv_vcf
        String chromosome
        String max_size
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu: 1,
                                   mem_gb: 8,
                                   disk_gb: 30,
                                   boot_disk_gb: 20,
                                   preemptible: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File filtered_vcf = "~{prefix}.cnv.~{chromosome}.vcf.gz"
        File filtered_vcf_index = "~{prefix}.cnv.~{chromosome}.vcf.gz.tbi"
    }

    command <<<
        set -euo pipefail
        echo "Subset to SVLEN <= max_cnv_size "
        bcftools view ~{gatk_sv_vcf} \
            -i 'INFO/SVLEN<=~{max_size}' \
            -o ~{prefix}.cnv.~{chromosome}.vcf.gz

        tabix -p vcf ~{prefix}.cnv.~{chromosome}.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_pipeline_docker
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

task subsetGATKSV {
	input {
        File gatk_sv_vcf
        File gatk_sv_vcf_idx
        File sample_list
        Int min_cnv_size
        Int max_cnv_size
        Int? max_ac
        String prefix
        String chromosome
        String scripts
        String array_validation_docker
        RuntimeAttr? runtime_attr_override
	}

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 32,
        disk_gb: 30,
        boot_disk_gb: 20,
        preemptible: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
        File subset_vcf = "~{prefix}.cnv.~{chromosome}.vcf.gz"
        File subset_vcf_index = "~{prefix}.cnv.~{chromosome}.vcf.gz.tbi"
	}

	command <<<
        set -euo pipefail
        echo "Subset to samples in sample list, contig of interest, DEL/DUP SVTYPEs, and SVLEN >= min_cnv_size "
        bcftools view ~{gatk_sv_vcf} \
            -r ~{chromosome} \
            -S ~{sample_list} \
            -i '(INFO/SVTYPE=="DEL" || INFO/SVTYPE=="DUP") && INFO/SVLEN>=~{min_cnv_size} && INFO/SVLEN<~{max_cnv_size}' \
            -O u \
            | bcftools view \
            --min-ac 1 \
            ~{"--max-ac " + max_ac} \
            -O z \
            -o ~{prefix}.cnv.~{chromosome}.vcf.gz

        tabix -p vcf ~{prefix}.cnv.~{chromosome}.vcf.gz
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

