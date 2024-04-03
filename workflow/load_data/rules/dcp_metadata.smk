"""
DCP Metadata
"""

if 'dcp_metadata' in config.keys():
    metadata_df = pd.read_table(config['dcp_metadata'],comment='#')
else:
    metadata_df = pd.DataFrame(columns=['study', 'filename'])

dcp_studies = set(metadata_df['study']).intersection(set(dataset_df['study']))


# rule download_dcp_tsv:
#     output:
#         tsv=out_dir / 'dcp_metadata' / '{study}' / 'dcp_metadata.tsv'
#     params:
#         url=lambda wildcards: metadata_df.query('study == @wildcards.study')['url'].values[0]
#     run:
#         shell("wget -O {output} '{params.url}'")


# rule download_dcp_all:
#     input:
#         expand(rules.download_dcp_tsv.output,study=dataset_df['study'])


rule add_dcp_metadata:
    input:
        zarr=rules.load_data_merge_study.output.zarr,
        dcp=lambda wildcards: metadata_df.query('study == @wildcards.study')['filename'].values[0],
        # dcp=rules.download_dcp_tsv.output.tsv,
    output:
        zarr=directory(out_dir / 'dcp_metadata' / '{study}.zarr'),
        stats=out_dir / 'dcp_metadata' / '{study}' / 'stats.tsv'
    params:
        id_cols=[
            'donor_organism.uuid',
            'donor_organism.biomaterial_core.biomaterial_id',
            'donor_organism.biomaterial_core.biomaterial_name',
            'cell_line.biomaterial_core.biomaterial_id',
            'organoid.biomaterial_core.biomaterial_id',
            'sample.biomaterial_core.biomaterial_id',
            'sequencing_input.biomaterial_core.biomaterial_id',
            'specimen_from_organism.uuid',
            'specimen_from_organism.biomaterial_core.biomaterial_id',
            'specimen_from_organism.biomaterial_core.biomaterial_name',
            'donor_id',
        ],
        metadata_cols=[
            'analysis_process',
            'biomaterial_core',
            'cell_suspension',
            'collection_protocol',
            'donor_organism',
            'enrichment_protocol',
            'library_preparation_protocol',
            'process',
            'protocol',
            'sequencing_protocol',
            'specimen_from_organism',
        ]
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs/')
    resources:
        mem_mb='10GB'
    script:
        '../scripts/add_dcp_metadata.py'


rule add_dcp_metadata_all:
    input:
        expand(rules.add_dcp_metadata.output,study=dcp_studies)
    localrule: True


def collect_stats(wildcards):
    return {
        study:
        expand(rules.add_dcp_metadata.output.stats,study=study)[0]
        for study in dcp_studies
    }


rule plot_stats:
    input:
        unpack(collect_stats)
    output:
        intersection_plot=image_dir / 'dcp_metadata' / 'id_intersection.png',
        intersection_stats=out_dir / 'dcp_metadata' / 'id_intersection.tsv',
    conda:
        get_env(config, 'plots', env_dir='../../../envs/')
    script:
        '../scripts/plot_stats.py'


rule dcp_metadata_all:
    input:
        rules.add_dcp_metadata_all.output,
        rules.plot_stats.output,
    localrule: True
