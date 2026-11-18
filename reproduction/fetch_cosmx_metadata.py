"""Fetch the pinned metadata member from Bruker's public ZIP."""

import argparse
from hashlib import sha256
import io
import json
from pathlib import Path
from zipfile import ZipFile

import requests

from source_asset import COSMX_ARCHIVE, COSMX_METADATA


class RemoteZip(io.RawIOBase):
    """Read only the ZIP ranges required to recover the metadata member."""

    def __init__(self):
        super().__init__()
        self.session = requests.Session()
        self.session.mount("https://", requests.adapters.HTTPAdapter(max_retries=3))
        with self.request(0, 0) as response:
            self.size = int(response.headers["Content-Range"].rsplit("/", 1)[1])
        self.position = 0

    def request(self, first, last):
        response = self.session.get(
            COSMX_ARCHIVE,
            headers={"Range": f"bytes={first}-{last}"},
            timeout=45,
            stream=True,
        )
        prefix = f"bytes {first}-{last}/"
        if response.status_code != 206 or not response.headers.get(
            "Content-Range", ""
        ).startswith(prefix):
            response.close()
            raise OSError(
                "Bruker ZIP does not support byte ranges. Use the manual download."
            )
        return response

    def tell(self):
        return self.position

    def seek(self, offset, whence=io.SEEK_SET):
        position = {
            io.SEEK_SET: 0,
            io.SEEK_CUR: self.position,
            io.SEEK_END: self.size,
        }[whence] + offset
        if position < 0:
            raise ValueError("Negative file position")
        self.position = position
        return position

    def seekable(self):
        return True

    def close(self):
        self.session.close()
        super().close()

    def read(self, size=-1):
        stop = self.size if size < 0 else min(self.size, self.position + size)
        if stop <= self.position:
            return b""
        with self.request(self.position, stop - 1) as response:
            data = response.content
        if len(data) != stop - self.position:
            raise OSError("Truncated response from the Bruker ZIP.")
        self.position = stop
        return data


def fetch(bundle, output):
    manifest = json.loads((bundle / "manifest.json").read_text())
    source = manifest["external_inputs"]["cosmx_metadata"]
    with RemoteZip() as remote, ZipFile(remote) as archive:
        (member,) = [n for n in archive.namelist() if Path(n).name == COSMX_METADATA]
        if archive.getinfo(member).file_size != source["bytes"]:
            raise ValueError("Unexpected CosMx metadata size.")
        data = archive.read(member)
    if sha256(data).hexdigest() != source["sha256"]:
        raise ValueError("Bruker's CosMx metadata does not match the frozen checksum.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream:
        stream.write(data)
    print(f"Verified CosMx metadata at {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", default=Path("."), type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    fetch(args.bundle, args.output)
