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
3. check_format: Asserts that the target file is readable as vcf.
"""

from pathlib import Path

from cyvcf2 import VCF

from ensembl.datacheck.functions.file_checks import file_exists, is_gz_text_file
from ensembl.datacheck.functions.io_utils import vcf_reader


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
