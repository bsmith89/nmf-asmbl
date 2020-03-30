MAX_THREADS = 30

noperiod='[^.]+'
integer='[0-9]*'
wildcard_constraints:
    radius=integer,
    key=noperiod,
    k=integer,


rule run_bcalm:
    output:
        fasta='data/core.bcalm-k{k}.fn',
        hdf5='data/core.bcalm-k{k}.h5',
    input:
        reads=[f'data/{stem}.m.proc.r1.fq.gz'
               for stem
               in ['SS01117', 'SS01120']]
    params:
        k=31,
        min_abund=2,
    threads: MAX_THREADS
    shell:
        """
        readslist=$(mktemp)
        tmpdir=$(mktemp -d)
        echo $tmpdir "(temporary directory for {output.fasta})" >&2

        for file in {input.reads}
        do
            realpath $file
        done > $readslist

        bcalm -nb-cores {threads} \
                -kmer-size {params.k} -abundance-min {params.min_abund} \
                -in $readslist -out-tmp $tmpdir -out $tmpdir/out
        echo "BCALM done; moving output files from $tmpdir"
        mv $tmpdir/out.unitigs.fa {output.fasta}
        mv $tmpdir/out.h5 {output.hdf5}
        """

rule convert_bcalm_to_gfa:
    output:
        gfa='data/{stem}.bcalm-k{k}.gfa'
    input:
        script='scripts/bcalm_to_gfa.py',
        fasta='data/{stem}.bcalm-k{k}.fn',
    shell:
        """
        {input.script} {input.fasta} {output.gfa} {wildcards.k}
        """

rule fetch_strain_graph_cloud:
    output: '{stem}.{key}.cloud-r{radius}.gfa'
    input:
        graph='{stem}.gfa',
        contig_list='{stem}.{key}.list',
    shell:
        """
        gfaview -S {wildcards.radius} -s {input.contig_list} {input.graph} > {output}
        """
