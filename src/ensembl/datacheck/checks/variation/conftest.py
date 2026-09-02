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

from ensembl.datacheck.checks.variation.types import Csq_subfield_spec

# Constants
DEFAULT_MAX_RANDOM_REGIONS = "1000"

# pytest fixtures
@pytest.fixture(scope="session")
def CSQ_SPECS() -> dict[str, Csq_subfield_spec]:
    """
    Dictionary of CSQ subfield specifications/requirements.

    Keys are CSQ subfield names, values are the subfield specs as `Csq_subfield_spec`.
    """
    return {
        "Allele": {"canbe_empty": False, "species": "all"},
        "Consequence": {"canbe_empty": False, "species": "all"},
        "Feature": {"species": "all"},
        "VARIANT_CLASS": {"canbe_empty": False, "species": "all"},
        "SPDI": {"canbe_empty": False, "species": "all"},
        "PUBMED": {"species": "all"},
        "VAR_SYNONYMS": {"species": "all"},
        "PHENOTYPES": {},
        "Conservation": {"species": "homo_sapiens"},
        "CADD_PHRED": {},
        "AA": {},
        "SIFT": {},
        "PolyPhen": {},
        "gnomAD_exomes_AF": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_afr": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_afr": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_afr": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_amr": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_amr": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_amr": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_asj": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_asj": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_asj": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_eas": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_eas": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_eas": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_fin": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_fin": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_fin": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_nfe": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_nfe": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_nfe": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_oth": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_oth": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_oth": {"species": "homo_sapiens"},
        "gnomAD_exomes_AF_sas": {"species": "homo_sapiens"},
        "gnomAD_exomes_AC_sas": {"species": "homo_sapiens"},
        "gnomAD_exomes_AN_sas": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_afr": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_afr": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_afr": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_amr": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_amr": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_amr": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_asj": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_asj": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_asj": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_eas": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_eas": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_eas": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_fin": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_fin": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_fin": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_nfe": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_nfe": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_nfe": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_oth": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_oth": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_oth": {"species": "homo_sapiens"},
        "gnomAD_genomes_AF_sas": {"species": "homo_sapiens"},
        "gnomAD_genomes_AC_sas": {"species": "homo_sapiens"},
        "gnomAD_genomes_AN_sas": {"species": "homo_sapiens"},
        "AF": {"species": "homo_sapiens"},
        "AFR_AF": {"species": "homo_sapiens"},
        "AMR_AF": {"species": "homo_sapiens"},
        "EAS_AF": {"species": "homo_sapiens"},
        "EUR_AF": {"species": "homo_sapiens"},
        "SAS_AF": {"species": "homo_sapiens"},
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
