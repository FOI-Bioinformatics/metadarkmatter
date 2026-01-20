"""
API clients for external services.

Provides clients for GTDB taxonomy database and other web APIs.
"""

from metadarkmatter.clients.gtdb import GTDBClient, GTDBGenome, GTDBQueryResult

__all__ = [
    "GTDBClient",
    "GTDBGenome",
    "GTDBQueryResult",
]
