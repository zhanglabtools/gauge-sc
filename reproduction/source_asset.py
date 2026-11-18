"""Shared panel inventory and checksum helpers for frozen results."""

import hashlib

PANEL_COUNTS = {2: 9, 3: 6, 4: 6, 5: 8, 6: 11, 7: 5}
SCHEMATICS = {"fig3a", "fig4a", "fig6a", "fig7a"}
COSMX_PANELS = {"fig6h", "fig6i"}
COSMX_METADATA = "Run5642_S3_Quarter_metadata_file.csv"
COSMX_ARCHIVE = "https://smi-public.objects.liquidweb.services/Quarter%20Brain.zip"
COSMX_ACCESS = (
    "https://brukerspatialbiology.com/resources/"
    "coronal-hippocampus-and-cortex-basic-data-files/"
)


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()
