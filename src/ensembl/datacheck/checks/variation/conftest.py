# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
variation/conftest.py

Pytest fixtures for variation-specific datachecks.

Fixtures provided:
1. variation_params: Merges variation default sampling parameters with any
   user-supplied --params values.
"""

import pytest

from ensembl.datacheck.checks.variation.types import Expectation_spec

# Constants
DEFAULT_MAX_RANDOM_REGIONS = "1000"

# pytest fixtures
@pytest.fixture(scope="session")
def CSQ_FIELD_EXPECTATIONS() -> dict[str, Expectation_spec]:
    return {
        "Allele": {"canbe_empty": False, "field_existance": "all"},
        "Consequence": {"canbe_empty": False, "field_existance": "all"},
        "Feature": {"field_existance": "all"},
        "VARIANT_CLASS": {"canbe_empty": False, "field_existance": "all"},
        "SPDI": {"canbe_empty": False, "field_existance": "all"},
        "PUBMED": {"field_existance": "all"},
        "VAR_SYNONYMS": {"field_existance": "all"},
        "PHENOTYPES": {},
        "Conservation": {"field_existance": "homo_sapiens"},
        "CADD_PHRED": {},
        "AA": {},
        "SIFT": {},
        "PolyPhen": {},
        "gnomAD_exomes_AF": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_afr": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_afr": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_afr": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_amr": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_amr": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_amr": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_asj": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_asj": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_asj": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_eas": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_eas": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_eas": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_fin": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_fin": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_fin": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_nfe": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_nfe": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_nfe": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_oth": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_oth": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_oth": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AF_sas": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AC_sas": {"field_existance": "homo_sapiens"},
        "gnomAD_exomes_AN_sas": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_afr": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_afr": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_afr": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_amr": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_amr": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_amr": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_asj": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_asj": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_asj": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_eas": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_eas": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_eas": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_fin": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_fin": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_fin": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_nfe": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_nfe": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_nfe": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_oth": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_oth": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_oth": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AF_sas": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AC_sas": {"field_existance": "homo_sapiens"},
        "gnomAD_genomes_AN_sas": {"field_existance": "homo_sapiens"},
        "AF": {"field_existance": "homo_sapiens"},
        "AFR_AF": {"field_existance": "homo_sapiens"},
        "AMR_AF": {"field_existance": "homo_sapiens"},
        "EAS_AF": {"field_existance": "homo_sapiens"},
        "EUR_AF": {"field_existance": "homo_sapiens"},
        "SAS_AF": {"field_existance": "homo_sapiens"},
    }


@pytest.fixture(scope="session")
def variation_params(params):
    """
    Variation-specific default parameters for source-file sampling checks.

    Args:
        params (dict): Parsed command-line params from the shared plugin.

    Returns:
        dict: User params merged over variation defaults.
    """
    resolved_params = {"max_random_regions": DEFAULT_MAX_RANDOM_REGIONS}
    resolved_params.update(params)
    return resolved_params
