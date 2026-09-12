#!/usr/bin/env python3
"""Generate the deterministic 64 MiB repeating/sparsely modified fixture."""

import argparse
import hashlib
from pathlib import Path


MIB = 1024 * 1024
COPIES = 64
PATCHES_PER_COPY = 32
SEED = b"SpeedSketch feasibility fixture v1"


def digest_bytes(label, size):
    return hashlib.shake_256(SEED + label).digest(size)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    base = digest_bytes(b"/base", MIB)
    sha256 = hashlib.sha256()
    with args.output.open("xb") as output:
        for copy in range(COPIES):
            block = bytearray(base)
            if copy:
                for patch in range(PATCHES_PER_COPY):
                    token = digest_bytes(copy.to_bytes(2, "little") + patch.to_bytes(2, "little"), 9)
                    offset = int.from_bytes(token[:8], "little") % MIB
                    block[offset] ^= token[8] or 1
            output.write(block)
            sha256.update(block)
    print(f"{sha256.hexdigest()}  {args.output}")


if __name__ == "__main__":
    main()
