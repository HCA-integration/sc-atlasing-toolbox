import scanpy as sc
import glob

outputs = glob.glob('test/out/integration/test/*.h5ad')
print(outputs)

for file in outputs:
    adata = sc.read(file)

    try:
        # Metadata
        assert 'dataset' in adata.uns
        assert adata.uns['dataset'] == 'test'
        assert 'methods' in adata.uns
        assert 'integration' in adata.uns
        for key in ['method', 'label_key', 'batch_key', 'output_type']:
            assert key in adata.uns['integration']

        # Unintegrated data
        assert 'counts' in adata.layers
        assert 'normcounts' in adata.layers

        # Output type specific outputs
        output_type = adata.uns['integration']['output_type']
        output_type = [output_type] if isinstance(output_type, str) else output_type

        if 'knn' in output_type:
            assert 'connectivities' in adata.obsp
            assert 'distances' in adata.obsp
        elif 'embed' in output_type:
            assert 'X_emb' in adata.obsm
            assert 'connectivities' in adata.obsp
            assert 'distances' in adata.obsp
        elif 'full' in output_type:
            assert 'corrected_counts' in adata.layers
            assert 'X_emb' in adata.obsm
            assert 'connectivities' in adata.obsp
            assert 'distances' in adata.obsp
        else:
            raise ValueError(f'Invalid output type {output_type}')

    except AssertionError as e:
        print(file)
        print(adata)
        raise e
