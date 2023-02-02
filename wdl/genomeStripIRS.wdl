version 1.0

import "Structs.wdl"

task genomeStripIRS {
    input {
        File input_file
        File genome
        File genome_index
        File genome_dict
        File array
        String prefix
        File gs_tarball
        String	array_validation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 150,
        disk_gb: 180,
        boot_disk_gb: 100,
        preemptible: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File vcf = "~{prefix}.irs.vcf.gz"
        File report = "~{prefix}.irs.report.dat.gz"
    }

    command <<<

        tar -xzvf ~{gs_tarball} -C /cromwell_root

        export SV_DIR=/cromwell_root/svtoolkit
        export classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"

        java -Xmx24g -cp $classpath \
        org.broadinstitute.sv.main.SVAnnotator \
        -A IntensityRankSum \
        -R ~{genome} \
        -vcf ~{input_file} \
        -O ~{prefix}.irs.vcf \
        -arrayIntensityFile ~{array} \
        -irsUseGenotypes true \
        -writeReport true \
        -reportFile ~{prefix}.irs.report.dat

        gzip ~{prefix}.irs.vcf
        gzip ~{prefix}.irs.report.dat
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