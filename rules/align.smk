from pathlib import Path

samples_dict = dict()
with open(expand_path(working_dir, "samples/samples.tsv")) as file_handle:
    for sample_line in file_handle.readlines():
        tab_sep_line = sample_line.strip().split("\t")
        assert len(tab_sep_line) == 2
        samples_dict[tab_sep_line[0]] = tab_sep_line[1]

samples = list(samples_dict.keys())

rule star_align:
    input:
        expand("{root_dir}samples/{samples}/{samples}.fastq.gz", root_dir=config["working_dir"], samples=samples)
    output:
        bam = expand_path(align_dir, "{sample}")
    shell:
        "echo 'dsa' > {output}"
