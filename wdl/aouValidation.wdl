version 1.0

import "genomeStripIRS.wdl" as gsirs

workflow ukbbArrayValidation {

    input {
        File samples_list
        String array_path
        File ids_corresp
        Array[String] contigs
        File genome
        File genome_index
        File genome_dict
        String gs_path
        String array_validation_docker
        RuntimeAttr? runtime_attr_override
    }

    Array[String] samples = transpose(read_tsv(samples_list))[0]

    scatter(sample in samples){
        #vcf_path = "gs://prod-drc-broad/Array_vcfs/hg38/"
        File sample_vcf = "~{array_path}/~{sample}.*.sorted.vcf.gz"
        File sample_idx = "~{sample_vcf}.tbi"

        call calculateLRR{
            input:
                input_vcf=sample_vcf,
                input_idx=sample_idx,
                prefix=sample,
                array_validation_docker=array_validation_docker,
                runtime_attr_override = runtime_attr_override
        }

    call mergeLRR{
        input:
            files=select_all(calculateLRR.array_lrr),
            prefix=prefix,
            array_validation_docker=array_validation_docker,
            runtime_attr_override = runtime_attr_override
    }

    scatter(contig in contigs){

        call subsetGATKSV{
            input:
                gatk_sv_vcf=gatk_sv_vcf,
                prefix=sample,
                chromosome=contig
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
                prefix=prefix,
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
        String prefix
        String array_validation_docker
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
        File lrr = "~{prefix}.lrr.gz"
        File array_lrr = "~{prefix}.lrr.exp"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r gs://alba-gatk-sv/array-validation/scripts .

        bcftools reheader --samples ~{samples_list} -o ~{prefix}.reheader.vcf.gz ~{input_vcf}

        bcftools query -H -f "%ID\t%CHROM\t%POS[\t%LRR]\n" ~{prefix}.reheader.vcf.gz | \
            gzip > ~{prefix}.lrr.gz

        echo "Calculating LRRs"
        python3 scripts/calculateLRR.py --input ~{prefix}.lrr.gz \
            --output ~{prefix}.lrr.exp
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
        echo "Copying scripts"
        gsutil -m cp -r gs://alba-gatk-sv/array-validation/scripts .

        echo "Subset only samples in sample list"
        bcftools view ~{gatk_sv_vcf} ~{chromosome} -S ~{sample_list} -O z -o gatk-sv.~{chomosome}.vcf.gz

        echo "Subset only DEL/DUP"
        bcftools view gatk-sv.~{chomosome}.vcf.gz | \
            grep -E "^#|DEL|DUP" | bcftools view -O z -o gatk-sv.cnv.~{chomosome}.vcf.gz
        tabix -p vcf gatk-sv.cnv.~{chomosome}.vcf.gz
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
        String array_validation_docker
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
        File lrr_files = "LRRfiles.fof"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r gs://alba-gatk-sv/array-validation/scripts .

        echo "Merging LRR files"
        echo "~{sep=" " files}" > LRRfiles.fof
        sed -i "s/ /\n/g" LRRfiles.fof
        Rscript scripts/mergeFiles.R -f LRRfiles.fof -o ~{prefix}.merged.report.dat
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