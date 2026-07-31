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

from ensembl.production.metadata.api.factories.genomes import GenomeFactory
from ensembl.production.metadata.api.models import Dataset, GenomeDataset, DatasetStatus, EnsemblRelease, ReleaseStatus, \
    Organism, Genome, OrganismGroup, OrganismGroupMember, Assembly, AssemblySequence, DatasetSource, GenomeRelease, \
    DatasetType
import sys
from typing import Any, Iterable, List, Optional

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing_extensions import override


class EnsemblDatacheckWarning(UserWarning):
    """
    Custom warning class for Ensembl data checks.

    Attributes:
        message (str): The warning message.
        file_name (str): The name of the file where the warning originated.
        function_name (str): The name of the function where the warning originated.
    """

    def __init__(self: Self, message: str, file_name: str, function_name: str) -> None:
        """
        Initializes the EnsemblDatacheckWarning with the given message, file name, and function name.

        Args:
            message (str): The warning message.
            file_name (str): The name of the file where the warning originated.
            function_name (str): The name of the function where the warning originated.
        """
        self.message: str = message
        self.file_name: str = file_name
        self.function_name: str = function_name

    @override
    def __str__(self: Self) -> str:
        """
        Returns a formatted string representation of the warning.

        Returns:
            str: Formatted warning message.
        """
        return f"Warning::{self.file_name}::{self.function_name}: {self.message}"


def get_genomes_from_metadata_db(
    db_url: str, release_name: Optional[str] = None,
    genome_uuids: Optional[List[str]] = None
) -> Iterable[dict[str, Any]]:
    """
    Fetch genome UUIDs from database based on release_name or explicit UUID list.

    Args:
        db_url: SQLAlchemy DB URL
        release_name: comma-separated list of release names
        genome_uuids: explicit genome UUIDs

    Returns:
        Iterable of dictionaries containing genome metadata including uuid
    """
    release_name_list: List[str] | None = release_name.split(',') if release_name else None
    genomes_iter: Iterable[dict[str, Any]] = GenomeFactory().get_genomes(
                metadata_db_uri=db_url,
                genome_uuid=genome_uuids,
                dataset_type="genebuild",
                dataset_names=["genebuild"],
                release_name=release_name_list,
                batch_size=0,  # Fetch all at once since we need to group them in memory
                columns=[
                    GenomeRelease.release_id.label("release_id"),
                    EnsemblRelease.name.label("release_name"),
                    EnsemblRelease.version.label("release_version"),
                    EnsemblRelease.status.label("release_status"),
                    EnsemblRelease.label.label("release_label"),
                    EnsemblRelease.is_current.label("release_is_current"),
                    Genome.genome_uuid.label("genome_uuid"),
                    Genome.production_name.label("species"),
                    Genome.created.label("genome_submitted"),
                    Genome.genebuild_date.label("genebuild_date"),
                    Organism.scientific_name.label("scientific_name"),
                    DatasetSource.name.label("db_name"),
                    Assembly.accession.label("assembly_accession"),
                    Assembly.name.label("assembly_name"),
                    Dataset.dataset_uuid.label("dataset_uuid"),
                    Dataset.status.label("dataset_status"),
                    Dataset.name.label("dataset_name"),
                    DatasetType.name.label("dataset_type"),
                    EnsemblRelease.name.label("genome_release"),
                    GenomeDataset.release_id.label("genome_dataset_release_id"),
                    GenomeDataset.is_current.label("genome_dataset_is_current"),
                ]
            )
    return genomes_iter
