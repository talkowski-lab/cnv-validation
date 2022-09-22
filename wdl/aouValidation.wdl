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
        String? min_cnv_size  # lower bound on SVLEN to evaluate against arrays. Default: 50000 (50kb)
        Int max_ac  # maximum allele count to evaluate against arrays

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

    Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]

    scatter (i in range(length(samples))) {
        call calculateLRR {
            input:
                input_vcf=array_vcfs[i],
                input_idx=array_vcfs[i] + ".tbi",
                contigs=contigs,
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
            runtime_attr_override = runtime_attr_merge_lrr
    }

    scatter (contig in contigs) {
        call subsetGATKSV {
            input:
                gatk_sv_vcf=gatk_sv_vcf,
                gatk_sv_vcf_idx="~{gatk_sv_vcf}.tbi",
                sample_list=write_lines(samples),
                prefix=prefix,
                max_ac=max_ac,
                min_cnv_size=select_first([min_cnv_size, "50000"]),
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
        Array[String] contigs
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
        set -euo pipefail
        echo "Copying scripts"
        gsutil -m cp -r ~{scripts} .

        echo "Reheader VCF and limit to primary contigs"
        bcftools view -r ~{sep="," contigs} -O u ~{input_vcf} \
            | bcftools reheader --samples ~{ids_corresp} \
            | bgzip -c > ~{sample}.reheader.vcf.gz

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
        String min_cnv_size
        Int max_ac
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
            --max-ac ~{max_ac}
            -i '(INFO/SVTYPE=="DEL" || INFO/SVTYPE=="DUP") && INFO/SVLEN>=~{min_cnv_size}' \
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

task mergeLRR {
	input {
        Array[File] files
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
        File merged_lrr = "~{prefix}.merged.lrr.exp.tsv"
	}

	command <<<
        set -euo pipefail

        echo "Merging LRR files"

        python3 <<CODE
        import pandas as pd
        lrr = ["~{sep='"," ' files}"]
        all_lrr = None
        for file in lrr:
            tmp = pd.read_table(file)
            if all_lrr is None:
                all_lrr = tmp.iloc[:,:4]
            all_lrr[tmp.columns[4]] = tmp.iloc[:,4]
        total_probes = len(all_lrr)
        all_lrr = all_lrr[np.all(all_lrr != ".", axis=1)]  # drop probe if any sample is missing LRR
        probes_full_data = len(all_lrr)
        dropped_probes = total_probes - probes_full_data
        print(f"Total probes: {total_probes}. Probes after dropping missing data: {probes_full_data}. Dropped probes: {dropped_probes}.")
        all_lrr.to_csv("~{prefix}.merged.lrr.exp.tsv", sep='\t', index=False, header=True)
        CODE
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