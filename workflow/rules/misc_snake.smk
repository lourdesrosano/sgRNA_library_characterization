# Helper functions required by the pipeline

import os.path
import sys
import pandas as pd
from snakemake.utils import validate


#########################################
###  Check config file for missing values
#########################################

fail_instantly = False


# Define error object and message
class Error(object):
    def __init__(self, key, name):
        self.__key = key
        self.__name = name

    def __add__(self, other):
        return self

    def __call__(self, wildcards=None):
        sys.exit(
            """
            ===============================================
            You have not specified '{}' for '{}'
            ===============================================
            """.format(
                self.__key, self.__name
            )
        )

    def __getitem__(self, value):
        return Error(key=self.__key, name=self.__name)


# Define config object
class Config(object):
    def __init__(self, kwargs, name="Config"):
        self.__name = name
        self.__members = {}
        for key, value in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=key)
            else:
                self.__members[key] = value

    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            if fail_instantly:
                sys.exit(
                    """
                    ===============================================
                    You have not specified '{}' for '{}'
                    ===============================================
                    """.format(
                        key, self.__name
                    )
                )
            else:
                return Error(key=key, name=self.__name)


# Check with the above class definitions if the config file contains all necessary values
config = Config(config)


###########################
###  READ INPUT FILES
###########################

## Import library map

# First, sanity check if was provided in the corresponding config
if isinstance(config["input_output"]["library_map"], Error):
    raise ValueError("'library_map' was not provided in the config!")

# Assumes header and that library IDs are provided in column "library_id" of the library map
library_map = pd.read_table(config["input_output"]["library_map"], header=0)
library_map = library_map.set_index("library_id", drop = False)

# Validate sample map using corresponding schema
validate(library_map, "../schema/library_map.schema.yaml")
# Column "library_id" is expected to contain the unique identifiers
library_ids = library_map["library_id"].tolist()


## Import TCGA samples

# First, sanity check if was provided in the corresponding config
if isinstance(config["input_output"]["tcga_samples"], Error):
    raise ValueError("'tcga_samples' was not provided in the config!")

# Assumes header and that TCGA sample IDs are provided in column "sample" of the file
tcga_samples = pd.read_table(config["input_output"]["tcga_samples"], header=0)
tcga_samples = tcga_samples.set_index("sample", drop = False)

# Validate sample file using corresponding schema
validate(tcga_samples, "../schema/tcga_samples.schema.yaml")
# Column "sample" is expected to contain the unique TCGA identifiers
tcga_sample_ids = tcga_samples["sample"].tolist()
