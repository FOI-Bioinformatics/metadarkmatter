"""
GTDB API client for querying genome accessions.

Provides access to the Genome Taxonomy Database (GTDB) API for
retrieving representative genome accessions by taxonomic group.
"""

from __future__ import annotations

import logging
import re
import time
from collections import Counter
from dataclasses import dataclass
from typing import Any, Self
from urllib.parse import quote

import httpx

logger = logging.getLogger(__name__)

from metadarkmatter.core.exceptions import MetadarkmatterError

GTDB_API_BASE = "https://gtdb-api.ecogenomic.org"

# Retry configuration defaults
DEFAULT_MAX_RETRIES = 3
DEFAULT_RETRY_DELAY = 1.0  # seconds
DEFAULT_RETRY_BACKOFF = 2.0  # exponential multiplier


class GTDBAPIError(MetadarkmatterError):
    """Error communicating with GTDB API."""

    def __init__(self, message: str, status_code: int | None = None):
        self.status_code = status_code
        suggestion = "Check your internet connection and try again."
        if status_code == 404:
            suggestion = "The taxon may not exist in GTDB. Check the spelling and format."
        elif status_code == 429:
            suggestion = "Rate limited. Wait a moment and try again."
        elif status_code and status_code >= 500:
            suggestion = "GTDB server error. Try again later."

        super().__init__(message=message, suggestion=suggestion)


class InvalidTaxonFormatError(MetadarkmatterError):
    """Invalid GTDB taxon format."""

    def __init__(self, taxon: str):
        self.taxon = taxon
        super().__init__(
            message=f"Invalid taxon format: '{taxon}'",
            suggestion=(
                "GTDB taxon format uses prefixes:\n"
                "  d__ = domain (e.g., d__Bacteria)\n"
                "  p__ = phylum (e.g., p__Proteobacteria)\n"
                "  c__ = class (e.g., c__Gammaproteobacteria)\n"
                "  o__ = order (e.g., o__Enterobacterales)\n"
                "  f__ = family (e.g., f__Enterobacteriaceae)\n"
                "  g__ = genus (e.g., g__Escherichia)\n"
                "  s__ = species (e.g., s__Escherichia coli)"
            ),
        )


@dataclass(frozen=True)
class GTDBGenome:
    """Single genome entry from GTDB.

    Attributes:
        accession: NCBI genome accession (e.g., GCF_000005845.2)
        gtdb_taxonomy: Full GTDB taxonomy string
        species: Species name
        genome_size: Genome size in base pairs (if available)
    """

    accession: str
    gtdb_taxonomy: str
    species: str
    genome_size: int | None = None


@dataclass(frozen=True)
class GTDBQueryResult:
    """Result of a GTDB taxonomy query.

    Attributes:
        taxon: Query taxon string
        genomes: Tuple of genome entries
        total_count: Total number of genomes
        genus_counts: Count of genomes per genus
        species_counts: Count of genomes per species
    """

    taxon: str
    genomes: tuple[GTDBGenome, ...]
    total_count: int
    genus_counts: dict[str, int]
    species_counts: dict[str, int]


class GTDBClient:
    """Client for GTDB API requests.

    Provides methods for querying genome information from the
    Genome Taxonomy Database (GTDB) API.

    Attributes:
        timeout: Request timeout in seconds
    """

    # Regex pattern for valid GTDB taxon format
    TAXON_PATTERN = re.compile(r"^[dpcofgs]__[\w\s]+$")

    def __init__(
        self,
        timeout: float = 60.0,
        max_retries: int = DEFAULT_MAX_RETRIES,
        retry_delay: float = DEFAULT_RETRY_DELAY,
        retry_backoff: float = DEFAULT_RETRY_BACKOFF,
    ):
        """Initialize GTDB client.

        Args:
            timeout: Request timeout in seconds.
            max_retries: Maximum retry attempts for transient failures.
            retry_delay: Initial delay between retries in seconds.
            retry_backoff: Exponential backoff multiplier for retries.
        """
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.retry_backoff = retry_backoff
        self._client: httpx.Client | None = None

    def _get_client(self) -> httpx.Client:
        """Get or create HTTP client."""
        if self._client is None:
            self._client = httpx.Client(
                base_url=GTDB_API_BASE,
                timeout=self.timeout,
                headers={"Accept": "application/json"},
            )
        return self._client

    def close(self) -> None:
        """Close the HTTP client."""
        if self._client is not None:
            self._client.close()
            self._client = None

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def __del__(self) -> None:
        """Ensure HTTP client is closed on garbage collection."""
        self.close()

    @classmethod
    def validate_taxon_format(cls, taxon: str) -> bool:
        """Check if taxon matches GTDB format.

        Args:
            taxon: Taxon string to validate

        Returns:
            True if valid GTDB format
        """
        return bool(cls.TAXON_PATTERN.match(taxon))

    def query_genomes(
        self,
        taxon: str,
        *,
        representatives_only: bool = True,
    ) -> GTDBQueryResult:
        """Query genomes for a taxonomic group.

        Args:
            taxon: GTDB taxon string (e.g., f__Enterobacteriaceae)
            representatives_only: Only return representative genomes

        Returns:
            GTDBQueryResult with genome list and statistics

        Raises:
            InvalidTaxonFormatError: If taxon format is invalid
            GTDBAPIError: If API request fails
        """
        if not self.validate_taxon_format(taxon):
            raise InvalidTaxonFormatError(taxon)

        # Use genomes-detail endpoint to get full taxonomy information
        # URL-encode taxon to handle special characters safely
        endpoint = f"/taxon/{quote(taxon, safe='')}/genomes-detail"
        params = {"sp_reps_only": str(representatives_only).lower()}

        response_data = self._get(endpoint, params)
        genomes = self._parse_genomes_detail(response_data)

        return GTDBQueryResult(
            taxon=taxon,
            genomes=genomes,
            total_count=len(genomes),
            genus_counts=self._count_by_rank(genomes, "genus"),
            species_counts=self._count_by_rank(genomes, "species"),
        )

    def _get(self, endpoint: str, params: dict[str, str]) -> Any:
        """Make GET request to GTDB API with retry logic.

        Implements exponential backoff for transient failures (5xx errors,
        connection errors, rate limiting).

        Args:
            endpoint: API endpoint path
            params: Query parameters

        Returns:
            JSON response data

        Raises:
            GTDBAPIError: If request fails after all retries
        """
        client = self._get_client()
        last_exception: Exception | None = None
        delay = self.retry_delay

        for attempt in range(self.max_retries + 1):
            try:
                response = client.get(endpoint, params=params)
                response.raise_for_status()
                return response.json()

            except httpx.HTTPStatusError as e:
                last_exception = e
                status_code = e.response.status_code

                # Don't retry client errors (4xx) except rate limiting (429)
                if 400 <= status_code < 500 and status_code != 429:
                    raise GTDBAPIError(
                        f"GTDB API request failed: {status_code}",
                        status_code=status_code,
                    ) from e

                # Handle rate limiting with Retry-After header
                if status_code == 429:
                    retry_after = e.response.headers.get("Retry-After")
                    if retry_after:
                        delay = float(retry_after)

                # Retry on server errors (5xx) or rate limiting (429)
                if attempt < self.max_retries:
                    logger.warning(
                        "GTDB API request failed (attempt %d/%d): %s. Retrying in %.1fs...",
                        attempt + 1,
                        self.max_retries + 1,
                        status_code,
                        delay,
                    )
                    time.sleep(delay)
                    delay *= self.retry_backoff

            except httpx.RequestError as e:
                last_exception = e
                # Retry on connection errors
                if attempt < self.max_retries:
                    logger.warning(
                        "GTDB API connection error (attempt %d/%d): %s. Retrying in %.1fs...",
                        attempt + 1,
                        self.max_retries + 1,
                        str(e),
                        delay,
                    )
                    time.sleep(delay)
                    delay *= self.retry_backoff

        # All retries exhausted
        if isinstance(last_exception, httpx.HTTPStatusError):
            raise GTDBAPIError(
                f"GTDB API request failed after {self.max_retries + 1} attempts: "
                f"{last_exception.response.status_code}",
                status_code=last_exception.response.status_code,
            ) from last_exception
        raise GTDBAPIError(
            f"GTDB API request failed after {self.max_retries + 1} attempts: "
            f"{last_exception}"
        ) from last_exception

    def _parse_genomes_detail(
        self, response_data: dict[str, Any] | list[Any]
    ) -> tuple[GTDBGenome, ...]:
        """Parse genome entries from genomes-detail API response.

        The genomes-detail endpoint returns {"rows": [...]} with each row
        containing individual taxonomy rank fields.

        Args:
            response_data: JSON response from GTDB API

        Returns:
            Tuple of GTDBGenome objects
        """
        # Handle the {"rows": [...]} format from genomes-detail endpoint
        if isinstance(response_data, dict):
            rows = response_data.get("rows", [])
        elif isinstance(response_data, list):
            rows = response_data
        else:
            return ()

        genomes = []
        for entry in rows:
            if not isinstance(entry, dict):
                continue

            # Extract accession from 'gid' field
            accession = entry.get("gid", "")

            # Build full taxonomy string from individual rank fields
            taxonomy_parts = []
            for rank in ["gtdbDomain", "gtdbPhylum", "gtdbClass", "gtdbOrder",
                         "gtdbFamily", "gtdbGenus", "gtdbSpecies"]:
                value = entry.get(rank, "")
                if value:
                    taxonomy_parts.append(value)

            gtdb_taxonomy = ";".join(taxonomy_parts)

            # Extract species name (remove s__ prefix)
            species_field = entry.get("gtdbSpecies", "")
            species = species_field.removeprefix("s__")

            # Genome size not available in genomes-detail response
            genome_size = entry.get("genomeSize")
            if genome_size is not None:
                try:
                    genome_size = int(genome_size)
                except (ValueError, TypeError):
                    genome_size = None

            if accession:
                genomes.append(
                    GTDBGenome(
                        accession=accession,
                        gtdb_taxonomy=gtdb_taxonomy,
                        species=species,
                        genome_size=genome_size,
                    )
                )

        return tuple(genomes)

    @staticmethod
    def _extract_species_from_taxonomy(taxonomy: str) -> str:
        """Extract species name from GTDB taxonomy string.

        Args:
            taxonomy: Full GTDB taxonomy string

        Returns:
            Species name or empty string
        """
        for part in taxonomy.split(";"):
            part = part.strip()
            if part.startswith("s__"):
                return part[3:]
        return ""

    @staticmethod
    def _extract_genus_from_taxonomy(taxonomy: str) -> str:
        """Extract genus name from GTDB taxonomy string.

        Args:
            taxonomy: Full GTDB taxonomy string

        Returns:
            Genus name or empty string
        """
        for part in taxonomy.split(";"):
            part = part.strip()
            if part.startswith("g__"):
                return part[3:]
        return ""

    def _count_by_rank(
        self, genomes: tuple[GTDBGenome, ...], rank: str
    ) -> dict[str, int]:
        """Count genomes by taxonomic rank.

        Args:
            genomes: Tuple of genome entries
            rank: Taxonomic rank ('genus' or 'species')

        Returns:
            Dictionary of name -> count
        """
        if rank == "genus":
            names = [
                self._extract_genus_from_taxonomy(g.gtdb_taxonomy) for g in genomes
            ]
        elif rank == "species":
            names = [g.species for g in genomes]
        else:
            names = []

        return dict(Counter(n for n in names if n))
