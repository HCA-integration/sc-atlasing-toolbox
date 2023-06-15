import pandas as pd


def load_dataset_df(file):
    df = pd.read_table(file,comment='#')
    df['url'] = df['url'].astype(str)
    df['subset'] = df['subset'].astype(str).str.split(',')
    df = df.explode('subset')
    return df