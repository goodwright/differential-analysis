process GEN_CONTRASTS {
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' : 'ubuntu:20.04' }"

    input:
    val input

    output:
    path "contrasts.csv", emit: file

    when:
    task.ext.when == null || task.ext.when

    script:

    output_string = "id,variable,reference,target,blocking\n"
    for(int i=0; i<input.size(); i++) {
        output_string += input[i][1] + "_" + input[i][2] + "_" + input[i][3] + "_" + input[i][4] + "," + input[i][1]  + "," + input[i][2]  + "," + input[i][3]  + "," + input[i][4] + "\n"
    }

    """
    echo "${output_string}" > contrasts.csv
    """
}
