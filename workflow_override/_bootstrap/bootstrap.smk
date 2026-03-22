import sys
sys.path.insert(0, "/home/ubuntu/scAtlasTb/workflow")
from pprint import pprint
from pathlib import Path
import pandas as pd
from snakemake.utils import min_version

from utils.pipeline import update_module_configs, config_for_module, update_input_files_per_dataset, update_file_for_module_param
from utils.ModuleConfig import ModuleConfig
from integration.IntegrationConfig import IntegrationConfig
from sample_representation.SampleRepresentationConfig import SampleRepresentationConfig
from preprocessing.PreprocessingConfig import PreprocessingConfig


min_version("6.0")
container: "docker://condaforge/mambaforge:latest"

configfile: "configs/outputs.yaml"
configfile: "configs/load_data/config.yaml"
configfile: "configs/exploration/config.yaml"

# params = pd.read_table('configs/modules.tsv',comment='#')
# params['submodules'] = params['submodules'].str.split(',')
# config = update_module_configs(config, params)

config_kwargs = {
    'batch_analysis': dict(
        config_params=['covariates', 'permute_covariates', 'n_permutations', 'sample'],
    ),
    'preprocessing': dict(
        config_params=['assemble', 'extra_hvgs'],
    ),
    'sample_representation': dict(
        parameters="/home/ubuntu/scAtlasTb/workflow/sample_representation/params.tsv",
        config_params=['methods', 'use_rep', 'var_mask', 'sample_key', 'norm_counts', 'raw_counts'],
        wildcard_names=['method', 'input_type', 'use_rep', 'var_mask'],
        rename_config_params={'methods': 'method'},
        explode_by=['method', 'use_rep', 'var_mask'],
    ),
    'split_data': dict(
        config_params=['key', 'values'],
        wildcard_names=['key', 'value'],
        rename_config_params={'values': 'value'},
        explode_by=['value'],
    ),
    'integration': dict(
        parameters="/home/ubuntu/scAtlasTb/workflow/integration/params.tsv",
        config_params=['methods', 'batch', 'label', 'var_mask'],
        wildcard_names=['method', 'batch', 'label', 'var_mask', 'output_type'],
        rename_config_params={'methods': 'method'},
        explode_by=['method', 'batch', 'label', 'var_mask'],
    ),
    'merge': dict(
        mandatory_wildcards=['dataset'],
    ),
    'collect': dict(
        mandatory_wildcards=['dataset'],
    ),
    'uncollect': dict(
        config_params=[
            'sep',
            'new_file_ids',
        ],
        wildcard_names=['new_file_id'],
        rename_config_params={'new_file_ids': 'new_file_id'},
        explode_by=['new_file_id'],
    ),
}

config_classes = {
    'integration': IntegrationConfig,
    'sample_representation': SampleRepresentationConfig,
    'preprocessing': PreprocessingConfig,
}

config['DATASETS'] = config.get('DATASETS', {})
default_datasets = config.get('defaults', {}).get('datasets', config['DATASETS'].keys())
for dataset, dataset_config in config['DATASETS'].items():
    for module_name in list(dataset_config.get('input', {}).keys()):
        config = update_input_files_per_dataset(
            dataset=dataset,
            module_name=module_name,
            config=config,
            config_class_map=config_classes,
            config_kwargs=config_kwargs,
        )


# TODO move to data loader and exploration
if 'dataset_meta' in config:
    config['dataset_meta'] = Path(config['dataset_meta']).resolve()

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

envvars:
    'HDF5_USE_FILE_LOCKING',
    'KMP_DUPLICATE_LIB_OK',

# Import modules
