import glob
import numpy as np
import scanpy as sc

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
        output_types = adata.uns['integration']['output_type']
        output_types = [output_types] if isinstance(output_types, str) else output_types

        if adata.uns['integration']['method'] == 'scanorama':
            assert len(output_types) == 2

        for ot in output_types:
            if ot not in ['knn', 'embed', 'full']:
                raise ValueError(f'Invalid output type {ot}')

        if 'knn' in output_types:
            assert 'connectivities' in adata.obsp
            assert 'distances' in adata.obsp
        if 'embed' in output_types:
            assert 'X_emb' in adata.obsm
        if 'full' in output_types:
            # check that counts are different from input
            adata_raw = adata.raw.to_adata()
            if adata.uns['integration']['method'] != 'unintegrated':
                assert adata.X.data.shape != adata_raw.X.data.shape or np.any(adata.X.data != adata_raw.X.data)
            assert 'X_pca' in adata.obsm

    except AssertionError as e:
        print(file)
        print(adata)
        raise e
