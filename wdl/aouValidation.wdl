version 1.0

import "genomeStripIRS.wdl" as gsirs
import "Structs.wdl"

workflow aouArrayValidation {

    input {
        File samples_list
        String array_path
        File ids_corresp
        Array[String] contigs
        File genome
        File genome_index
        File genome_dict
        File gatk_sv_vcf
        String gs_path
        String array_validation_docker
        RuntimeAttr? runtime_attr_override
    }

    Array[String] samples = transpose(read_tsv(samples_list))[0]

    scatter(sample in samples){
        File sample_vcf = "~{array_path}/~{sample}.*.sorted.vcf.gz"
        File sample_idx = "~{sample_vcf}.tbi"

        call calculateLRR{
            input:
                input_vcf=sample_vcf,
                input_idx=sample_idx,
                sample=sample,
                array_validation_docker=array_validation_docker,
                runtime_attr_override = runtime_attr_override
        }
    }
    call mergeLRR{
        input:
            files=select_all(calculateLRR.array_lrr),
            array_validation_docker=array_validation_docker,
            runtime_attr_override = runtime_attr_override
    }

    scatter(contig in contigs){

        call subsetGATKSV{
            input:
                gatk_sv_vcf=gatk_sv_vcf,
                gatk_sv_vcf_idx="~{gatk_sv_vcf}.tbi",
                chromosome=contig,
                array_validation_docker=array_validation_docker,
                runtime_attr_override = runtime_attr_override
        }

        call gsirs.genomeStripIRS as genomeStripIRS{
            input:
                input_file=subsetGATKSV.merged_vcf,
                genome=genome,
                genome_index=genome_index,
                genome_dict=genome_dict,
                array=mergeLRR.merged_lrr,
                samples_list=samples_list,
                gs_path=gs_path,
                array_validation_docker=array_validation_docker,
                runtime_attr_override = runtime_attr_override
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
        File samples_list
        String sample
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
        echo "Reheader VCF"
        bcftools reheader --samples ~{samples_list} -o ~{sample}.reheader.vcf.gz ~{input_vcf}

        bcftools query -H -f "%ID\t%CHROM\t%POS[\t%LRR]\n" ~{sample}.reheader.vcf.gz | \
            gzip > ~{sample}.lrr.gz

        echo "Calculating LRRs"
        python3 calculateLRR.py --input ~{sample}.lrr.gz \
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
        String chromosome
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
        File merged_vcf = "merged_gatk-sv.vcf.gz"
	}

	command <<<
        echo "Subset only samples in sample list"
        bcftools view ~{gatk_sv_vcf} ~{chromosome} -S ~{sample_list} -O z -o gatk-sv.~{chromosome}.vcf.gz

        echo "Subset only DEL/DUP"
        bcftools view gatk-sv.~{chromosome}.vcf.gz | \
            grep -E "^#|DEL|DUP" | bcftools view -O z -o gatk-sv.cnv.~{chromosome}.vcf.gz
        tabix -p vcf gatk-sv.cnv.~{chromosome}.vcf.gz
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
        File merged_lrr = "aou.merged.report.dat"
        File lrr_files = "LRRfiles.fof"
	}

	command <<<
        echo "Merging LRR files"
        echo "~{sep=" " files}" > LRRfiles.fof
        sed -i "s/ /\n/g" LRRfiles.fof
        Rscript mergeFiles.R -f LRRfiles.fof -o aou.merged.report.dat
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