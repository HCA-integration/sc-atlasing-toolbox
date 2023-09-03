from .misc import ifelse


def get_env(
    config,
    env_name,
    mode=None,
    env_dir='../envs',
    gpu_env=None,
):
    if mode is None:
        mode = config['env_mode'] if 'env_mode' in config else 'local'
    if gpu_env is None:
        gpu_env = env_name
    gpu_env = ifelse(
        'use_gpu' not in config or config['use_gpu'],
        _if=env_name,
        _else=gpu_env
    )

    if mode == 'from_yaml':
        return f'{env_dir}/{env_name}.yaml'
    elif mode == 'local':
        return env_name
    else:
        raise ValueError(f'Unknown environment mode: {mode}')

