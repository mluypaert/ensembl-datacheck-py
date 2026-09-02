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

from ensembl.datacheck.checks.variation.types import Csq_subfield_spec
from ensembl.datacheck.checks.vcf import *
from ensembl.datacheck.functions.vcf_utils import (
    get_vcf_variant_count_by_chr,
    get_vcf_variant_count,
    parse_CSQ_format,
    vcf_reader
)

# module level marker to select test per filetype and dataset type
pytestmark = [
    pytest.mark.dataset_type('short_variants'),
    pytest.mark.file_extension("vcf.gz"),
    pytest.mark.file_extension("vcf")
]


# Header checks
def check_csq_in_header(subtests: pytest.Subtests, target_file: Path, csq_specs_species_filtered: dict[str, Csq_subfield_spec]):
    """
    Check that the VCF INFO header contains the CSQ field.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the VCF info header does not contain the CSQ field.
    """

    reader = vcf_reader(target_file)

    assert reader.get_header_type('CSQ'), "CSQ field not found in the VCF (INFO) header."

    # Parse the CSQ format described in the field description
    csq_fields = parse_CSQ_format(target_file)

    for field,spec in csq_specs_species_filtered.items():
        with subtests.test(f"Evaluating presence of CSQ subfield '{field}'", field=field, spec=spec):
            assert field in csq_fields, f"CSQ field '{field}' not found in the CSQ fromat string."


def check_source_in_header(target_file: Path):
    """
    Check that the VCF header contains the source field.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the VCF header does not contain the source field.
    """

    reader = vcf_reader(target_file)

    assert reader.get_header_type('source'), "source field not found in the VCF header."


# Source comparison checks
def check_variant_count_source_comparison(target_file: Path, source_file: Path | None):
    """
    Compare target VCF variant count with source VCF variant count.

    The check asserts:
    - source_file is provided
    - at least 90% of the source variant count is found in the target
    - source variant count does not exceed target variant count

    Args:
        target_file: Path to target VCF file.
        source_file: Path to source VCF file (required).

    Raises:
        `AssertionError` if:
          * inputs are missing/invalid
          * count ratio is below threshold
          * target count exceeds source count
    """
    assert source_file is not None, "A source file is required (--source-file)."

    target_variant_count: int | None = get_vcf_variant_count(target_file)
    source_variant_count: int | None = get_vcf_variant_count(source_file)

    assert target_variant_count is not None
    assert source_variant_count is not None
    assert (target_variant_count / source_variant_count) > 0.90, "Target file variant count is less than 90% of the source file variant count."
    assert (target_variant_count / source_variant_count) <= 1, "Target file contains more variants than the source file."


def check_per_chr_variant_count_source_comparison(subtests: pytest.Subtests, target_file: Path, source_file: Path | None):
    """
    Compare target VCF variant counts with source VCF variant counts per chromosome.
    Only comparses common chromosomes between target and source.

    The check asserts:
    - source_file is provided
    - for each of the common chromosomes between target and source:
        * at least 95% of the source variant count is found in the target
        * source variant count does not exceed target variant count

    Args:
        target_file: Path to target VCF file.
        source_file: Path to source VCF file (required).

    Raises:
        `AssertionError` if:
          * inputs are missing/invalid
          * Count ratio is below threshold.
          * Target count exceeds source count for any chr.
    """
    assert source_file is not None, "A source file is required (--source-file)."

    target_chrs: list[str]
    with vcf_reader(target_file) as reader:
        target_chrs = reader.seqnames
    source_chrs: list[str]
    with vcf_reader(source_file) as reader:
        source_chrs = reader.seqnames

    common_chrs = list(set(target_chrs) & set(source_chrs))

    target_variant_counts = get_vcf_variant_count_by_chr(str(target_file))
    source_variant_counts = get_vcf_variant_count_by_chr(str(source_file))

    assert target_variant_counts is not None
    assert source_variant_counts is not None

    for chr in common_chrs:
        # Note: Not all chr will be present in source VCF (vcf_prepper rename some of them),
        #       nor in the target vcf (vcf_prepper remove some chr variants).
        if chr in target_variant_counts and chr in source_variant_counts:
            with subtests.test(f"Comparing variant counts for chr '{chr}'", chr=chr):
                recovery_ratio = target_variant_counts[chr] / source_variant_counts[chr]
                assert recovery_ratio > 0.90, f"ADVISORY: Target file variant count is only {recovery_ratio * 100}%  (<90%) of the source file variant count for chr {chr}."
                assert recovery_ratio <= 1, f"Target file contains more variants than the source file for chr {chr}."


# Content checks
def check_subsample_csq_content(csq_specs_species_filtered: dict[str, Csq_subfield_spec], subtests: pytest.Subtests, target_variants_subsample: dict):
    """
    Check that the CSQ field is populated as expected
    for a random subset of variants from the target file.

    Args:
        target_variants_subsample (dict): A subsample of variants from the target file to check.
        csq_specs_species_filtered: Dictionary of CSQ subfield specifications, filtered to fields relevant for the input species only.

    Raises:
        AssertionError: If the CSQ field is not populated as expected in any of the tested variants.
    """
    assert target_variants_subsample is not None and len(target_variants_subsample) > 0, "Failed to sample variants from target file."

    for field,spec in csq_specs_species_filtered.items():
        with subtests.test(f"Evaluating CSQ subfield '{field}' content", field=field, spec=spec):
            canbe_empty = spec.get('canbe_empty', True)

            csq_field_cnt = 0
            for variant_id in target_variants_subsample:
                first_csq = target_variants_subsample[variant_id]['csqs'][0]
                if first_csq.get(field, "") != "":
                    csq_field_cnt += 1

            subsample_size = len(target_variants_subsample)
            if not canbe_empty:
                assert csq_field_cnt == subsample_size, f"Required CSQ field '{field}' is missing in {subsample_size - csq_field_cnt} out of {subsample_size} sampled variants."
            else:
                assert csq_field_cnt > 0, f"ADVISORY: Optional CSQ field '{field}' not found in any of the sampled variants."
