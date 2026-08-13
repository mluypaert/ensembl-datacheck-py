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

from gzip import GzipFile
import os
from pathlib import Path


def file_exists(target_file):
    """
    Check if the file exists at the given path.

    Args:
        target_file (str): The path to the file.

    Returns:
        bool: True if the file exists, False otherwise.
    """
    if not target_file:
        return False
    return os.path.exists(target_file)


def file_size(target_file):
    """
    Get the size of the file at the given path.

    Args:
        target_file (str): The path to the file.

    Returns:
        int or None: The size of the file in bytes,
                     or None if the file does not exist.
    """
    if not target_file or not os.path.exists(target_file):
        return None
    return os.path.getsize(target_file)


def is_text_file(target_file: str | Path) -> bool:
    """
    Check if the file at the given path is a text file.

    Args:
        target_file: The path to the file.

    Returns:
        bool: True if the file is a text file, False otherwise.
    """
    try:
        with open(target_file, 'r', encoding='utf-8') as file:
            file.read()
        return True
    except Exception:
        return False


def is_gz_text_file(target_file: str | Path) -> bool:
    """
    Check if the file at the given path is a gzipped text file.

    Args:
        target_file: The path to the file.

    Returns:
        bool: True if the file is a text file, False otherwise.
    """
    try:
        with open(target_file, 'rb') as file:
            decompressed_file = GzipFile(fileobj=file)
            decompressed_file.read()
        return True
    except Exception:
        return False
