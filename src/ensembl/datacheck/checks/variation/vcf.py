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
variation/vcf.py

This module performs variation-specific vcf checks.
"""

from pathlib import Path
import pytest

from ensembl.datacheck.checks.vcf import *  # noqa: W0401
from ensembl.datacheck.functions.io_utils import vcf_reader

# module level marker to select test per filetype and dataset type
pytestmark = [
    pytest.mark.dataset_type('short_variants'),
    pytest.mark.file_extension("vcf.gz"),
    pytest.mark.file_extension("vcf")
]


def check_csq_in_header(target_file: str | Path):
    """
    Check that the VCF INFO header contains the CSQ field.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the VCF info header does not contain the CSQ field.
    """

    reader = vcf_reader(target_file)

    assert reader.get_header_type('CSQ'), "CSQ field not found in the VCF (INFO) header."


def check_source_in_header(target_file: str | Path):
    """
    Check that the VCF header contains the source field.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the VCF header does not contain the source field.
    """

    reader = vcf_reader(target_file)

    assert reader.get_header_type('source'), "source field not found in the VCF header."
