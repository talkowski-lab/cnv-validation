version 1.0

import "genomeStripIRS.wdl" as gsirs
import "Structs.wdl"

workflow aouArrayValidation {

    input {
        Array[String] samples
        Array[File] array_vcfs
        File gatk_sv_vcf
        File ids_corresp
        String prefix

        File primary_contigs_fai
        File genome
        File genome_index
        File genome_dict

        File scripts
        String gs_path
        String array_validation_docker

        RuntimeAttr? runtime_attr_calculate_lrr
        RuntimeAttr? runtime_attr_merge_lrr
        RuntimeAttr? runtime_attr_subset_gatk_sv
        RuntimeAttr? runtime_attr_genome_strip_irs
    }

    scatter (i in range(length(samples))) {
        call calculateLRR {
            input:
                input_vcf=array_vcfs[i],
                input_idx=array_vcfs[i] + ".tbi",
                sample=samples[i],
                ids_corresp=ids_corresp,
                array_validation_docker=array_validation_docker,
                scripts=scripts,
                runtime_attr_override = runtime_attr_calculate_lrr
        }
    }
    call mergeLRR {
        input:
            files=select_all(calculateLRR.array_lrr),
            array_validation_docker=array_validation_docker,
            prefix=prefix,
            scripts=scripts,
            runtime_attr_override = runtime_attr_merge_lrr
    }

    Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]
    scatter (contig in contigs) {
        call subsetGATKSV {
            input:
                gatk_sv_vcf=gatk_sv_vcf,
                gatk_sv_vcf_idx="~{gatk_sv_vcf}.tbi",
                sample_list=write_lines(samples),
                prefix=prefix,
                chromosome=contig,
                array_validation_docker=array_validation_docker,
                scripts=scripts,
                runtime_attr_override = runtime_attr_subset_gatk_sv
        }

        call gsirs.genomeStripIRS as genomeStripIRS{
            input:
                input_file=subsetGATKSV.subset_vcf,
                prefix=prefix,
                genome=genome,
                genome_index=genome_index,
                genome_dict=genome_dict,
                array=mergeLRR.merged_lrr,
                samples_list=write_lines(samples),
                gs_path=gs_path,
                array_validation_docker=array_validation_docker,
                runtime_attr_override = runtime_attr_genome_strip_irs
            }
        }

    output {
        Array[File] irs_vcf = genomeStripIRS.vcf
        Array[File] irs_report = genomeStripIRS.report
    }

}

task calculateLRR {
	input {
        File input_vcf
        File input_idx
        File ids_corresp
        String sample
        String scripts
        String array_validation_docker
        RuntimeAttr? runtime_attr_override

	}

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 16,
        disk_gb: 30,
        boot_disk_gb: 12,
        preemptible: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
        File lrr = "~{sample}.lrr.gz"
        File array_lrr = "~{sample}.lrr.exp"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r ~{scripts} .

        echo "Reheader VCF"
        bcftools reheader --samples ~{ids_corresp} -o ~{sample}.reheader.vcf.gz ~{input_vcf}

        bcftools query -H -f "%ID\t%CHROM\t%POS[\t%LRR]\n" ~{sample}.reheader.vcf.gz | \
            gzip > ~{sample}.lrr.gz

        echo "Calculating LRRs"
        python3 scripts/calculateLRR.py --input ~{sample}.lrr.gz \
            --output ~{sample}.lrr.exp
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
        echo "Subset only samples in sample list"
        bcftools view ~{gatk_sv_vcf} ~{chromosome} -S ~{sample_list} -O z -o ~{prefix}.~{chromosome}.vcf.gz

        echo "Subset only DEL/DUP"
        bcftools view ~{prefix}.~{chromosome}.vcf.gz | \
            grep -E "^#|DEL|DUP" | bcftools view -O z -o ~{prefix}.cnv.~{chromosome}.vcf.gz
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

task mergeLRR {
	input {
        Array[File] files
        String prefix
        String scripts
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
        File merged_lrr = "~{prefix}.merged.report.dat"
        File lrr_files = "~{prefix}.LRRfiles.fof"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r ~{scripts} .

        echo "Merging LRR files"
        echo "~{sep=" " files}" > ~{prefix}.LRRfiles.fof
        sed -i "s/ /\n/g" ~{prefix}.LRRfiles.fof
        Rscript scripts/mergeFiles.R -f ~{prefix}.LRRfiles.fof -o ~{prefix}.merged.report.dat
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