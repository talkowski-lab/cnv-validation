version 1.0

import "genomeStripIRS.wdl" as gsirs

workflow ukbbArrayValidation {

    input {
#        File array_bcf
#        String prefix
        File array_fof
        String array_validation_docker
        File samples_list
        Array[String] contigs
        String gcnv_dir
        File header
        File genome
        File genome_index
        File genome_dict
        String gs_path
    }

    Array[String] array_files = transpose(read_tsv(array_fof))[0]

    scatter(array_bcf in array_files){

        String prefix = basename(array_bcf, ".bcf")

        call getVCF{
            input:
                array_bcf=array_bcf,
                samples_list=samples_list,
                prefix=prefix,
                array_validation_docker=array_validation_docker
        }

        scatter (contig in contigs) {
            call calculateLRR{
                input:
                    ori_vcf=getVCF.ori_vcf,
                    ori_vcf_idx=getVCF.ori_vcf_idx,
                    samples_list=getVCF.samples,
                    prefix=prefix,
                    array_validation_docker=array_validation_docker,
                    chromosome=contig
            }

            File gcnv_file = "~{gcnv_dir}/ukbb_199811_gCNV.~{contig}.bed.gz"

            call gCNVbed2vcf{
                input:
                    samples_list=getVCF.samples,
                    gcnv_file=gcnv_file,
                    header=header,
                    prefix=prefix,
                    array_validation_docker=array_validation_docker,
                    chromosome=contig
            }
        }

        call mergeLRR{
            input:
                files=select_all(calculateLRR.array_lrr),
                prefix=prefix,
                array_validation_docker=array_validation_docker
        }

        call mergeVCF{
            input:
                files=select_all(gCNVbed2vcf.gcnv_vcf),
                prefix=prefix,
                array_validation_docker=array_validation_docker
        }

        call gsirs.genomeStripIRS as genomeStripIRS{
            input:
                input_file=mergeVCF.merged_vcf,
                genome=genome,
                genome_index=genome_index,
                genome_dict=genome_dict,
                array=mergeLRR.merged_lrr,
                samples_list=getVCF.samples,
                prefix=prefix,
                gs_path=gs_path,
                array_validation_docker=array_validation_docker
        }
    }

    output {
            Array[File] irs_vcf = genomeStripIRS.vcf
            Array[File] irs_report = genomeStripIRS.report
    }

}

task getVCF {
	input {
        File array_bcf
        File samples_list
        String prefix
        String array_validation_docker
	}

	output {
        File ori_vcf = "~{prefix}.vcf.gz"
        File ori_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File samples = "samples.list"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r gs://alba-gatk-sv/array-validation/scripts .

        echo "Reformatting UKBB VCF file"
        bcftools view ~{array_bcf} -O z -o ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz

        grep -wf <(bcftools query -l ~{prefix}.vcf.gz) ~{samples_list} > samples.list
	>>>

	runtime {
        memory: "32 GiB"
        disks: "local-disk 80 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: array_validation_docker
        bootDiskSizeGb: 40
    }
}

task calculateLRR {
	input {
        File ori_vcf
        File ori_vcf_idx
        File samples_list
        String prefix
        String chromosome
        String array_validation_docker
	}

	output {
        File lrr = "~{prefix}.~{chromosome}.lrr.gz"
        File array_lrr = "~{prefix}.~{chromosome}.lrr.exp"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r gs://alba-gatk-sv/array-validation/scripts .

        bcftools view ~{ori_vcf} ~{chromosome} -S ~{samples_list} --force-samples | \
            bcftools query -H -f "%ID\t%CHROM\t%POS[\t%LRR]\n" | \
        gzip > ~{prefix}.~{chromosome}.lrr.gz

        echo "Calculating LRRs"
        python3 scripts/calculateLRR.py --input ~{prefix}.~{chromosome}.lrr.gz \
            --output ~{prefix}.~{chromosome}.lrr.exp \
            --samples ~{samples_list}
	>>>

	runtime {
        memory: "32 GiB"
        disks: "local-disk 50 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: array_validation_docker
        bootDiskSizeGb: 30
    }
}

task mergeLRR {
	input {
        Array[File] files
        String prefix
        String array_validation_docker
	}

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
        memory: "64 GiB"
        disks: "local-disk 80 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: array_validation_docker
        bootDiskSizeGb: 30
    }
}


task gCNVbed2vcf {
	input {
        File samples_list
        File gcnv_file
        File header
        String prefix
        String chromosome
        String array_validation_docker
	}

	output {
        File gcnv_vcf = "~{prefix}.~{chromosome}.gCNV.vcf"
	}

	command <<<
        echo "Copying scripts"
        gsutil -m cp -r gs://alba-gatk-sv/array-validation/scripts .

        echo "Subsetting gCNV file"
        Rscript scripts/cnvSubset.R -s ~{samples_list} -i ~{gcnv_file} -o ~{prefix}_~{chromosome}.bed

        echo "Reformatting gCNV BED to VCF"
        python scripts/bed2vcf.py --bed ~{prefix}_~{chromosome}.bed \
            --samples ~{samples_list} \
            --header ~{header} \
            --vcf ~{prefix}.~{chromosome}.gCNV.vcf

	>>>

	runtime {
        memory: "32 GiB"
        disks: "local-disk 30 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: array_validation_docker
        bootDiskSizeGb: 20
    }
}

task mergeVCF {
	input {
        Array[File] files
        String prefix
        String array_validation_docker
	}

	output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
	}

	command <<<
        echo "Concatenating VCF files"
        bcftools concat ~{sep=" " files} | bcftools sort -O z -o ~{prefix}.merged.vcf.gz
	>>>

	runtime {
        memory: "32 GiB"
        disks: "local-disk 30 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: array_validation_docker
        bootDiskSizeGb: 20
    }
}
