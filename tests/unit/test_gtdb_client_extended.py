"""Extended tests for GTDB API client.

Tests retry logic, HTTP error handling, response parsing edge cases,
and context manager behavior using mocked httpx calls.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import httpx
import pytest

from metadarkmatter.clients.gtdb import (
    GTDBAPIError,
    GTDBClient,
    GTDBGenome,
    InvalidTaxonFormatError,
)


# =============================================================================
# GTDBAPIError tests
# =============================================================================


class TestGTDBAPIError:
    """Tests for GTDBAPIError with various status codes."""

    def test_error_without_status_code(self) -> None:
        """Should use default suggestion when no status code is provided."""
        error = GTDBAPIError("Connection failed")
        assert error.status_code is None
        assert "internet connection" in str(error)

    def test_error_with_404_status(self) -> None:
        """Should suggest checking taxon name for 404 errors."""
        error = GTDBAPIError("Not found", status_code=404)
        assert error.status_code == 404
        assert "may not exist" in str(error)

    def test_error_with_429_status(self) -> None:
        """Should suggest waiting for rate limit errors."""
        error = GTDBAPIError("Rate limited", status_code=429)
        assert error.status_code == 429
        assert "Rate limited" in str(error)

    def test_error_with_500_status(self) -> None:
        """Should suggest trying later for server errors."""
        error = GTDBAPIError("Server error", status_code=500)
        assert error.status_code == 500
        assert "server error" in str(error).lower()

    def test_error_with_502_status(self) -> None:
        """Should suggest trying later for 502 gateway errors."""
        error = GTDBAPIError("Bad gateway", status_code=502)
        assert error.status_code == 502
        assert "server error" in str(error).lower()

    def test_error_with_400_status(self) -> None:
        """Should use default suggestion for generic 400 errors."""
        error = GTDBAPIError("Bad request", status_code=400)
        assert error.status_code == 400
        assert "internet connection" in str(error)


# =============================================================================
# GTDBClient._get() retry logic tests
# =============================================================================


class TestGTDBClientRetryLogic:
    """Tests for the HTTP GET retry mechanism."""

    def _make_http_status_error(
        self, status_code: int, headers: dict[str, str] | None = None
    ) -> httpx.HTTPStatusError:
        """Helper to create a mock HTTPStatusError."""
        mock_response = MagicMock(spec=httpx.Response)
        mock_response.status_code = status_code
        mock_response.headers = headers or {}
        mock_request = MagicMock(spec=httpx.Request)
        return httpx.HTTPStatusError(
            f"HTTP {status_code}",
            request=mock_request,
            response=mock_response,
        )

    def test_successful_request(self) -> None:
        """Should return JSON on first successful attempt."""
        client = GTDBClient(max_retries=2, retry_delay=0.01)
        mock_response = MagicMock()
        mock_response.json.return_value = {"rows": []}
        mock_response.raise_for_status = MagicMock()

        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.return_value = mock_response
        client._client = mock_http_client

        result = client._get("/test", {})
        assert result == {"rows": []}
        assert mock_http_client.get.call_count == 1

    def test_retry_on_500_error(self) -> None:
        """Should retry on 500 server error and succeed on second attempt."""
        client = GTDBClient(max_retries=2, retry_delay=0.01)

        error_500 = self._make_http_status_error(500)
        success_response = MagicMock()
        success_response.json.return_value = {"rows": [{"gid": "GCF_001"}]}
        success_response.raise_for_status = MagicMock()

        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = [error_500, success_response]
        client._client = mock_http_client

        with patch("metadarkmatter.clients.gtdb.time.sleep"):
            result = client._get("/test", {})

        assert result == {"rows": [{"gid": "GCF_001"}]}
        assert mock_http_client.get.call_count == 2

    def test_no_retry_on_404_error(self) -> None:
        """Should not retry on 404 client error (raises immediately)."""
        client = GTDBClient(max_retries=3, retry_delay=0.01)

        error_404 = self._make_http_status_error(404)
        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = error_404
        client._client = mock_http_client

        with pytest.raises(GTDBAPIError) as exc_info:
            client._get("/test", {})

        assert exc_info.value.status_code == 404
        assert mock_http_client.get.call_count == 1

    def test_no_retry_on_400_error(self) -> None:
        """Should not retry on 400 bad request error."""
        client = GTDBClient(max_retries=3, retry_delay=0.01)

        error_400 = self._make_http_status_error(400)
        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = error_400
        client._client = mock_http_client

        with pytest.raises(GTDBAPIError) as exc_info:
            client._get("/test", {})

        assert exc_info.value.status_code == 400
        assert mock_http_client.get.call_count == 1

    def test_retry_on_429_rate_limit(self) -> None:
        """Should retry on 429 rate limit with Retry-After header."""
        client = GTDBClient(max_retries=2, retry_delay=0.01)

        error_429 = self._make_http_status_error(
            429, headers={"Retry-After": "0.5"}
        )
        success_response = MagicMock()
        success_response.json.return_value = {"rows": []}
        success_response.raise_for_status = MagicMock()

        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = [error_429, success_response]
        client._client = mock_http_client

        with patch("metadarkmatter.clients.gtdb.time.sleep") as mock_sleep:
            result = client._get("/test", {})

        assert result == {"rows": []}
        # Should have used the Retry-After value as delay
        mock_sleep.assert_called_once()
        assert mock_sleep.call_args[0][0] == 0.5

    def test_retry_on_connection_error(self) -> None:
        """Should retry on connection errors."""
        client = GTDBClient(max_retries=2, retry_delay=0.01)

        request_error = httpx.RequestError("Connection refused")
        success_response = MagicMock()
        success_response.json.return_value = {"rows": []}
        success_response.raise_for_status = MagicMock()

        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = [request_error, success_response]
        client._client = mock_http_client

        with patch("metadarkmatter.clients.gtdb.time.sleep"):
            result = client._get("/test", {})

        assert result == {"rows": []}
        assert mock_http_client.get.call_count == 2

    def test_all_retries_exhausted_http_error(self) -> None:
        """Should raise GTDBAPIError after all retries exhausted (HTTP error)."""
        client = GTDBClient(max_retries=1, retry_delay=0.01)

        error_500 = self._make_http_status_error(500)
        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = error_500
        client._client = mock_http_client

        with patch("metadarkmatter.clients.gtdb.time.sleep"):
            with pytest.raises(GTDBAPIError, match="after 2 attempts"):
                client._get("/test", {})

        # 1 initial + 1 retry = 2 calls
        assert mock_http_client.get.call_count == 2

    def test_all_retries_exhausted_connection_error(self) -> None:
        """Should raise GTDBAPIError after all retries exhausted (connection error)."""
        client = GTDBClient(max_retries=1, retry_delay=0.01)

        request_error = httpx.RequestError("Connection timeout")
        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = request_error
        client._client = mock_http_client

        with patch("metadarkmatter.clients.gtdb.time.sleep"):
            with pytest.raises(GTDBAPIError, match="after 2 attempts"):
                client._get("/test", {})

    def test_exponential_backoff(self) -> None:
        """Should apply exponential backoff between retries."""
        client = GTDBClient(
            max_retries=3, retry_delay=1.0, retry_backoff=2.0
        )

        error_500 = self._make_http_status_error(500)
        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.side_effect = error_500
        client._client = mock_http_client

        with patch("metadarkmatter.clients.gtdb.time.sleep") as mock_sleep:
            with pytest.raises(GTDBAPIError):
                client._get("/test", {})

        # Delays: 1.0, 2.0, 4.0 (3 retries with backoff=2.0)
        delays = [call[0][0] for call in mock_sleep.call_args_list]
        assert len(delays) == 3
        assert delays[0] == pytest.approx(1.0)
        assert delays[1] == pytest.approx(2.0)
        assert delays[2] == pytest.approx(4.0)


# =============================================================================
# GTDBClient._parse_genomes_detail() tests
# =============================================================================


class TestParseGenomesDetail:
    """Tests for parsing the genomes-detail API response."""

    def test_parse_dict_with_rows(self) -> None:
        """Should parse response with rows key."""
        client = GTDBClient()
        response = {
            "rows": [
                {
                    "gid": "GCF_000005845.2",
                    "gtdbDomain": "d__Bacteria",
                    "gtdbPhylum": "p__Pseudomonadota",
                    "gtdbClass": "c__Gammaproteobacteria",
                    "gtdbOrder": "o__Enterobacterales",
                    "gtdbFamily": "f__Enterobacteriaceae",
                    "gtdbGenus": "g__Escherichia",
                    "gtdbSpecies": "s__Escherichia coli",
                    "genomeSize": 4641652,
                },
            ]
        }

        genomes = client._parse_genomes_detail(response)

        assert len(genomes) == 1
        assert genomes[0].accession == "GCF_000005845.2"
        assert genomes[0].species == "Escherichia coli"
        assert genomes[0].genome_size == 4641652
        assert "d__Bacteria" in genomes[0].gtdb_taxonomy
        assert "s__Escherichia coli" in genomes[0].gtdb_taxonomy

    def test_parse_list_response(self) -> None:
        """Should parse response that is a plain list of entries."""
        client = GTDBClient()
        response = [
            {
                "gid": "GCF_001",
                "gtdbSpecies": "s__Test species",
            },
        ]

        genomes = client._parse_genomes_detail(response)

        assert len(genomes) == 1
        assert genomes[0].accession == "GCF_001"
        assert genomes[0].species == "Test species"

    def test_parse_empty_rows(self) -> None:
        """Should return empty tuple for empty rows."""
        client = GTDBClient()
        result = client._parse_genomes_detail({"rows": []})
        assert result == ()

    def test_parse_non_dict_entry_skipped(self) -> None:
        """Should skip non-dict entries in the rows list."""
        client = GTDBClient()
        response = {"rows": ["not_a_dict", {"gid": "GCF_001", "gtdbSpecies": "s__X"}]}

        genomes = client._parse_genomes_detail(response)

        assert len(genomes) == 1
        assert genomes[0].accession == "GCF_001"

    def test_parse_entry_without_accession_skipped(self) -> None:
        """Should skip entries with empty or missing gid."""
        client = GTDBClient()
        response = {
            "rows": [
                {"gid": "", "gtdbSpecies": "s__Species A"},
                {"gid": "GCF_002", "gtdbSpecies": "s__Species B"},
            ]
        }

        genomes = client._parse_genomes_detail(response)

        assert len(genomes) == 1
        assert genomes[0].accession == "GCF_002"

    def test_parse_invalid_genome_size(self) -> None:
        """Should set genome_size to None for non-numeric values."""
        client = GTDBClient()
        response = {
            "rows": [
                {
                    "gid": "GCF_003",
                    "gtdbSpecies": "s__Test",
                    "genomeSize": "not_a_number",
                },
            ]
        }

        genomes = client._parse_genomes_detail(response)

        assert len(genomes) == 1
        assert genomes[0].genome_size is None

    def test_parse_non_dict_non_list_response(self) -> None:
        """Should return empty tuple for unexpected response type."""
        client = GTDBClient()
        result = client._parse_genomes_detail("unexpected_string")  # type: ignore[arg-type]
        assert result == ()

    def test_parse_missing_taxonomy_fields(self) -> None:
        """Should handle entries with missing taxonomy rank fields."""
        client = GTDBClient()
        response = {
            "rows": [
                {
                    "gid": "GCF_004",
                    "gtdbDomain": "d__Bacteria",
                    # Missing most taxonomy fields
                },
            ]
        }

        genomes = client._parse_genomes_detail(response)

        assert len(genomes) == 1
        assert genomes[0].gtdb_taxonomy == "d__Bacteria"
        assert genomes[0].species == ""


# =============================================================================
# GTDBClient._count_by_rank() tests
# =============================================================================


class TestCountByRank:
    """Tests for counting genomes by taxonomic rank."""

    def test_count_by_genus(self) -> None:
        """Should count genomes per genus."""
        client = GTDBClient()
        genomes = (
            GTDBGenome(
                accession="GCF_001",
                gtdb_taxonomy="g__Escherichia;s__Escherichia coli",
                species="Escherichia coli",
            ),
            GTDBGenome(
                accession="GCF_002",
                gtdb_taxonomy="g__Escherichia;s__Escherichia fergusonii",
                species="Escherichia fergusonii",
            ),
            GTDBGenome(
                accession="GCF_003",
                gtdb_taxonomy="g__Salmonella;s__Salmonella enterica",
                species="Salmonella enterica",
            ),
        )

        counts = client._count_by_rank(genomes, "genus")

        assert counts == {"Escherichia": 2, "Salmonella": 1}

    def test_count_by_species(self) -> None:
        """Should count genomes per species."""
        client = GTDBClient()
        genomes = (
            GTDBGenome(accession="GCF_001", gtdb_taxonomy="", species="E. coli"),
            GTDBGenome(accession="GCF_002", gtdb_taxonomy="", species="E. coli"),
            GTDBGenome(accession="GCF_003", gtdb_taxonomy="", species="S. enterica"),
        )

        counts = client._count_by_rank(genomes, "species")

        assert counts == {"E. coli": 2, "S. enterica": 1}

    def test_count_unknown_rank(self) -> None:
        """Should return empty dict for unknown rank."""
        client = GTDBClient()
        genomes = (
            GTDBGenome(accession="GCF_001", gtdb_taxonomy="", species="X"),
        )

        counts = client._count_by_rank(genomes, "phylum")
        assert counts == {}

    def test_count_excludes_empty_names(self) -> None:
        """Should exclude genomes with empty names from counts."""
        client = GTDBClient()
        genomes = (
            GTDBGenome(accession="GCF_001", gtdb_taxonomy="", species="E. coli"),
            GTDBGenome(accession="GCF_002", gtdb_taxonomy="", species=""),
        )

        counts = client._count_by_rank(genomes, "species")
        assert counts == {"E. coli": 1}


# =============================================================================
# GTDBClient.query_genomes() integration-level tests (mocked HTTP)
# =============================================================================


class TestQueryGenomes:
    """Tests for the query_genomes() method with mocked HTTP responses."""

    def test_query_returns_result_with_counts(self) -> None:
        """Should return a GTDBQueryResult with genus and species counts."""
        client = GTDBClient(max_retries=0, retry_delay=0.01)

        mock_response = MagicMock()
        mock_response.json.return_value = {
            "rows": [
                {
                    "gid": "GCF_001",
                    "gtdbGenus": "g__Francisella",
                    "gtdbSpecies": "s__Francisella tularensis",
                },
                {
                    "gid": "GCF_002",
                    "gtdbGenus": "g__Francisella",
                    "gtdbSpecies": "s__Francisella novicida",
                },
            ]
        }
        mock_response.raise_for_status = MagicMock()

        mock_http_client = MagicMock(spec=httpx.Client)
        mock_http_client.get.return_value = mock_response
        client._client = mock_http_client

        result = client.query_genomes("f__Francisellaceae")

        assert result.taxon == "f__Francisellaceae"
        assert result.total_count == 2
        assert len(result.genomes) == 2
        assert result.genus_counts.get("Francisella") == 2
        assert result.species_counts.get("Francisella tularensis") == 1

    def test_query_invalid_taxon_raises(self) -> None:
        """Should raise InvalidTaxonFormatError for bad taxon format."""
        client = GTDBClient()

        with pytest.raises(InvalidTaxonFormatError):
            client.query_genomes("BadFormat")


# =============================================================================
# GTDBClient context manager and lifecycle tests
# =============================================================================


class TestGTDBClientLifecycle:
    """Tests for client lifecycle management."""

    def test_close_idempotent(self) -> None:
        """Should safely call close() multiple times."""
        client = GTDBClient()
        # Force client creation
        _ = client._get_client()
        assert client._client is not None

        client.close()
        assert client._client is None

        # Second close should not raise
        client.close()
        assert client._client is None

    def test_del_closes_client(self) -> None:
        """Should close client on garbage collection."""
        client = GTDBClient()
        _ = client._get_client()
        assert client._client is not None

        client.__del__()
        assert client._client is None

    def test_get_client_creates_once(self) -> None:
        """Should reuse existing client on repeated calls."""
        client = GTDBClient()
        first = client._get_client()
        second = client._get_client()
        assert first is second
        client.close()

    def test_custom_retry_parameters(self) -> None:
        """Should accept custom retry configuration."""
        client = GTDBClient(
            max_retries=5,
            retry_delay=2.0,
            retry_backoff=3.0,
        )
        assert client.max_retries == 5
        assert client.retry_delay == 2.0
        assert client.retry_backoff == 3.0
