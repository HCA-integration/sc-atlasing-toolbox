"""
DCP Metadata
"""

if 'dcp_metadata' in config.keys():
    metadata_df = pd.read_table(config['dcp_metadata'],comment='#')
else:
    metadata_df = pd.DataFrame(columns=['study', 'filename'])

studies = set(metadata_df['study']).intersection(set(dataset_df['study']))

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


rule save_obs:
    input:
        zarr=rules.filter.output.zarr
    output:
        obs=out_dir / 'dcp_metadata' / '{study}' / 'obs.tsv'
    conda:
        '../../../envs/scanpy.yaml'
    script:
        '../scripts/save_obs.py'


rule obs_merge_dcp:
    input:
        obs=rules.save_obs.output.obs,
        dcp=lambda wildcards: metadata_df.query('study == @wildcards.study')['filename'].values[0],
        # dcp=rules.download_dcp_tsv.output.tsv,
    output:
        obs=out_dir / 'dcp_metadata' / '{study}' / 'obs_merged.tsv',
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
            'sequencing_protocol.instrument_manufacturer_model',
            'sequencing_protocol.instrument_manufacturer_model.ontology',
            'sequencing_protocol.instrument_manufacturer_model.ontology_label',
            'sequencing_protocol.method',
            'sequencing_protocol.method.ontology',
            'sequencing_protocol.method.ontology_label',
            'library_preparation_protocol.cell_barcode.barcode_read',
            'library_preparation_protocol.cell_barcode.barcode_offset',
            'library_preparation_protocol.cell_barcode.barcode_length',
            'library_preparation_protocol.input_nucleic_acid_molecule.ontology',
            'library_preparation_protocol.input_nucleic_acid_molecule.ontology_label',
            'library_preparation_protocol.library_construction_method.ontology',
            'library_preparation_protocol.library_construction_method.ontology_label',
            'library_preparation_protocol.end_bias',
            'library_preparation_protocol.strand',
            'library_preparation_protocol.umi_barcode.barcode_read',
            'library_preparation_protocol.umi_barcode.barcode_offset',
            'library_preparation_protocol.umi_barcode.barcode_length',
            'collection_protocol.method.ontology',
            'collection_protocol.method.ontology_label',
            'enrichment_protocol.method.ontology',
            'enrichment_protocol.method.ontology_label',
            'cell_suspension.biomaterial_core.biomaterial_id'
        ]
    conda:
        '../../../envs/scanpy.yaml'
    script:
        '../scripts/obs_merge_dcp.py'


rule obs_merge_dcp_all:
    input:
        expand(rules.obs_merge_dcp.output,study=studies)


def collect_stats(wildcards):
    return {
        study:
        expand(rules.obs_merge_dcp.output.stats,study=study)[0]
        for study in studies
    }


rule plot_stats:
    input:
        unpack(collect_stats)
    output:
        intersection=out_dir / 'dcp_metadata' / 'plots' / 'intersection.png'
    conda:
        '../../../envs/plots.yaml'
    script:
        '../scripts/plot_stats.py'


rule dcp_metadata_all:
    input:
        rules.obs_merge_dcp_all.output,
        rules.plot_stats.output.intersection
