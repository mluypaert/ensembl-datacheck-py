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
This module performs checks on the new Ensembl metadata.

Checks performed:
1. Database connection check (Error): Ensures the database session is available.
2. Tables not empty check (Error): Ensures that all tables contain data.
3. Datasets attached to a release (Error): Ensures all released datasets are attached to a release.


These checks are run to verify the integrity and completeness of the metadata in the new Ensembl metadata database.
"""

import pytest
from typing import Any, List, Tuple

from ensembl.production.metadata.api.models import Dataset, GenomeDataset, DatasetStatus, EnsemblRelease, ReleaseStatus
from sqlalchemy import or_, func

from ensembl.datacheck.functions.db_checks import (
    database_connection_check, tables_not_empty_check
)


@pytest.mark.usefixtures("db_session")
def check_database(db_session: pytest.Session) -> None:
    """
    Check if the database connection is established.

    Args:
        db_session (sqlalchemy.orm.Session): The database session.

    Raises:
        AssertionError: If the database session is not available.
    """
    assert database_connection_check(db_session), "Database session is not available"

@pytest.mark.usefixtures("db_session")
def check_tables(db_session: pytest.Session) -> None:
    """
    Check that all tables in the database are not empty.

    Args:
        db_session (sqlalchemy.orm.Session): The database session.

    Raises:
        AssertionError: If any table is empty.
    """
    result, message = tables_not_empty_check(db_session)
    assert result, message


@pytest.mark.usefixtures("db_session")
def test_released_datasets_have_released_releases(db_session: pytest.Session) -> None:
    """
    Check that all Released datasets are properly attached to Released releases.

    Performs two checks:
    1. Released datasets must have a genome_dataset entry with a release_id
    2. Released datasets must only be attached to Released releases

    Args:
        db_session (sqlalchemy.orm.Session): The database session.

    Raises:
        AssertionError: If any Released dataset is not properly attached to a Released release.
    """
    # Check 1: Find Released datasets without a genome_dataset or with null release_id
    datasets_without_release: List[Any] = (
        db_session.query(Dataset)
        .outerjoin(Dataset.genome_datasets)
        .filter(
            Dataset.status == DatasetStatus.RELEASED,
            or_(
                GenomeDataset.dataset_id == None,  # No genome_dataset entry
                GenomeDataset.release_id == None  # Has genome_dataset but no release_id
            )
        )
        .all()
    )

    if datasets_without_release:
        no_release_dataset_ids: List[str] = [ds.dataset_uuid for ds in datasets_without_release]
        assert False, f"Found {len(datasets_without_release)} Released datasets without a release: {no_release_dataset_ids}"

    # Check 2: Find Released datasets attached to non-Released releases
    datasets_with_unreleased_release: List[Tuple[Dataset, GenomeDataset, EnsemblRelease]] = (
        db_session.query(Dataset, GenomeDataset, EnsemblRelease)
        .join(Dataset.genome_datasets)
        .join(GenomeDataset.ensembl_release)
        .filter(
            Dataset.status == DatasetStatus.RELEASED,
            EnsemblRelease.status != ReleaseStatus.RELEASED
        )
        .all()
    )

    if datasets_with_unreleased_release:
        unreleased_dataset_ids: List[str] = [ds[0].dataset_uuid for ds in datasets_with_unreleased_release]
        assert False, f"Found {len(datasets_with_unreleased_release)} Released datasets attached to non-Released releases: {unreleased_dataset_ids}"
