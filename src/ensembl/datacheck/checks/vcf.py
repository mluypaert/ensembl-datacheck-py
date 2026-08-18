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
vcf.py

This module performs generic vcf checks.

Checks performed:
1. check_exist: Asserts that the target vcf path exists.
2. check_is_gz_txt_file: Asserts that the target file is a gzipped text file.
3. check_has_index_file: Asserts that an accompanying index file can be found.
4. check_format: Asserts that the target file is readable as vcf.
"""

from pathlib import Path

from cyvcf2 import VCF
from pysam import HTSFile, TabixFile, VariantFile

from ensembl.datacheck.functions.file_checks import file_exists, is_gz_text_file
from ensembl.datacheck.functions.vcf_utils import is_bgzf_compressed_file, is_vcf_file, vcf_reader


def check_exist(target_file: Path):
    """
    Check that the target file exists on disk.

    Args:
        target_file: Path to the target file.

    Raises:
        AssertionError: If the target file is missing.
    """
    assert file_exists(target_file), "The target VCF file does not exist."


def check_is_gz_txt_file(target_file: str | Path):
    """
    Check that the file is a gzipped text file.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the file is not identified as a gzipped text file.
    """
    assert is_gz_text_file(target_file), "The file is not identified as a gzipped text file."


def check_is_bgzf_vcf_file(target_file: str | Path):
    """
    Check that the file is a bgzipped vcf file.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the file is not identified as a bgzipped vcf file.
    """
    try:
        file: HTSFile
        with VariantFile(str(target_file), mode="r") as file:
            # Assess compression
            assert (
                is_bgzf_compressed_file(file)
            ), f"File is not identified as a bgzipped file: {target_file}"

            # Assess format
            assert (
                is_vcf_file(file)
            ), f"File is not identified as a vcf file: {target_file}"

    except (OSError, ValueError, NotImplementedError) as ex:
        raise AssertionError(
            f"Exception caught during VariantFile assessment of {target_file}:"
            + f" {ex}."
        ) from ex


def check_has_valid_index_file(target_file: str | Path):
    """
    Check that an accompanying index file can be found and read appropriately.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the expected index file is not found or not readable
            as a TabixFile index.
    """

    csi_path = str(target_file) + '.csi'
    tbi_path = str(target_file) + '.tbi'

    index_file: str
    if file_exists(csi_path):
        index_file = csi_path
    elif file_exists(tbi_path):
        index_file = tbi_path
    else:
        raise AssertionError(
            f"Could not find index file for VCF file {target_file}."
        )

    try:
        with TabixFile(filename=str(target_file), index=index_file, mode="r"):
            pass
    except (ValueError) as ex:
        raise AssertionError(
            f"Missing index file: {index_file}"
        ) from ex
    except IOError as ex:
        raise AssertionError(
            f"Failed to open index file: {index_file}"
        ) from ex


def check_format(target_file: str | Path):
    """
    Check that the target file is readable as vcf.

    Args:
        target_file: Path to the target file.

    Raises:
        AssertionError: If the file is missing, unreadable, or not vcf.
    """
    assert file_exists(target_file), "The target file does not exist."

    reader = None
    try:
        reader = vcf_reader(target_file)
        assert reader is not None, "Could not open target file as VCF."
        assert isinstance(reader, VCF), "The target file is not recognised as VCF."
    except Exception as exc:
        raise AssertionError(
            f"Failed to read the target file as VCF: {exc}"
        ) from exc
    finally:
        if reader is not None:
            reader.close()


def check_header(target_file: str | Path):
    """
    Check that the VCF header has all fields required for bcftools processing.

    Args:
        target_file: The path to the file.

    Raises:
        AssertionError: If the file is not identified as a gzipped text file.
    """

    reader = vcf_reader(target_file)

    assert reader.get_header_type('fileformat') is not None, "The fileformat field is missing from the VCF header."

    column_header_line = str(reader.raw_header.rstrip('\n').split('\n')[-1])
    assert column_header_line.startswith('#'), "Column header line does not start with '#'."

    columns = column_header_line[1:].split('\t')

    assert columns[0] == 'CHROM', "Column header 'CHROM' not found at expected position."
    assert columns[1] == 'POS', "Column header 'POS' not found at expected position."
    assert columns[2] == 'ID', "Column header 'ID' not found at expected position."
    assert columns[3] == 'REF', "Column header 'REF' not found at expected position."
    assert columns[4] == 'ALT', "Column header 'ALT' not found at expected position."
