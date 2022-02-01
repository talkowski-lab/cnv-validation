version 1.0

task genomeStripIRS {
    input {
        File input_file
        File genome
        File genome_index
        File genome_dict
        File array
        File samples_list
        String prefix
        String gs_path
        String	array_validation_docker
        String chromosome
    }

    output {
        File vcf = "${prefix}.~{chromosome}.irs.vcf.gz"
        File report = "${prefix}.~{chromosome}.irs.report.dat.gz"
    }

    command <<<

        gsutil -m cp -r ~{gs_path} /cromwell_root/

        export SV_DIR=/cromwell_root/svtoolkit
        export classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"

        java -Xmx8g -cp $classpath \
        org.broadinstitute.sv.main.SVAnnotator \
        -A IntensityRankSum \
        -R ~{genome} \
        -vcf ~{input_file} \
        -O ~{prefix}.~{chromosome}.irs.vcf \
        -arrayIntensityFile ~{array} \
        -sample ~{samples_list} \
        -irsSampleTag SAMPLES \
        -writeReport true \
        -reportFile ~{prefix}.~{chromosome}.irs.report.dat

        gzip ~{prefix}.~{chromosome}.irs.vcf
        gzip ~{prefix}.~{chromosome}.irs.report.dat
	>>>

    runtime {
        memory: "24 GiB"
        disks: "local-disk 20 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: array_validation_docker
        bootDiskSizeGb: 20
    }
}
