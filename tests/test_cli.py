#!/usr/bin/env python3
"""Black-box checks for the SpeedSketch CLI; no third-party packages needed."""

import argparse
import hashlib
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BINARY = Path("./speedsketch").resolve()
TIMEOUT = 30.0
SCHEMES = ("ss-g", "ss-g-s", "ss-x", "od-g", "od-x")
OD_ARGS = ("--features", "12", "--super-features", "3")
SCHEMA_KEYS = {
    "schema_version",
    "input",
    "scheme",
    "pipeline",
    "verify",
    "features",
    "super_features",
    "input_bytes",
    "chunks",
    "unique_chunks",
    "duplicate_chunks",
    "similar_chunks",
    "delta_chunks",
    "self_compressed_chunks",
    "raw_chunks",
    "after_dedup_bytes",
    "after_delta_bytes",
    "estimated_final_payload_bytes",
    "source_read_bytes",
    "base_read_bytes",
    "estimated_drr",
    "dce_pct",
    "exclusion_rate_pct",
    "false_positive_rate_pct",
    "wall_seconds",
    "chunking_seconds",
    "dedup_seconds",
    "sketch_seconds",
    "encoding_seconds",
    "verification_seconds",
    "wall_mib_s",
    "chunking_mib_s",
    "dedup_mib_s",
    "sketch_mib_s",
    "encoding_mib_s",
    "verification_checks",
    "verification_failures",
}
COUNT_KEYS = {
    "input_bytes",
    "chunks",
    "unique_chunks",
    "duplicate_chunks",
    "similar_chunks",
    "delta_chunks",
    "self_compressed_chunks",
    "raw_chunks",
    "after_dedup_bytes",
    "after_delta_bytes",
    "estimated_final_payload_bytes",
    "source_read_bytes",
    "base_read_bytes",
    "verification_checks",
    "verification_failures",
}
RATE_KEYS = {
    "estimated_drr",
    "dce_pct",
    "exclusion_rate_pct",
    "false_positive_rate_pct",
}
TIMING_KEYS = {key for key in SCHEMA_KEYS if key.endswith("_seconds") or key.endswith("_mib_s")}


def fixed_random(size):
    return hashlib.shake_256(b"SpeedSketch CLI fixture v1").digest(size)


def repeating(size):
    word = b"SpeedSketch-repeating-fixture\0"
    return (word * ((size + len(word) - 1) // len(word)))[:size]


class SpeedSketchCliTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not BINARY.is_file():
            raise FileNotFoundError(f"binary not found: {BINARY}; pass --binary PATH")
        cls._tmp = tempfile.TemporaryDirectory(prefix="speedsketch-cli-")
        cls.tmp = Path(cls._tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def fixture(self, name, data):
        path = self.tmp / name
        path.write_bytes(data)
        return path

    def invoke(self, *args, expected=0):
        command = [str(BINARY), *map(str, args)]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=TIMEOUT)
        except subprocess.TimeoutExpired as error:
            self.fail(f"command timed out after {TIMEOUT}s: {' '.join(command)}\n{error}")
        self.assertEqual(
            expected,
            result.returncode,
            f"command: {' '.join(command)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def run_json(self, path, scheme="ss-g", pipeline="sequential", no_verify=False):
        args = ["--input", path, "--scheme", scheme, "--json"]
        if pipeline is not None:
            args += ["--pipeline", pipeline]
        if scheme.startswith("od-"):
            args += OD_ARGS
        if no_verify:
            args.append("--no-verify")
        result = self.invoke(*args)
        try:
            payload = json.loads(result.stdout)
        except json.JSONDecodeError as error:
            self.fail(f"stdout is not one JSON object: {error}\n{result.stdout}")
        self.assertEqual(SCHEMA_KEYS, set(payload), payload)
        self.assert_json_shape(payload)
        if scheme.startswith("od-"):
            self.assertEqual((12, 3), (payload["features"], payload["super_features"]))
        else:
            self.assertIsNone(payload["features"])
            self.assertIsNone(payload["super_features"])
        return payload

    def assert_json_shape(self, payload):
        self.assertIs(type(payload["schema_version"]), int)
        self.assertIs(type(payload["input"]), str)
        self.assertIs(type(payload["verify"]), bool)
        self.assertIn(payload["scheme"], SCHEMES)
        self.assertIn(payload["pipeline"], ("sequential", "parallel"))
        for key in COUNT_KEYS:
            self.assertIs(type(payload[key]), int, key)
            self.assertGreaterEqual(payload[key], 0, key)
        for key in RATE_KEYS | TIMING_KEYS:
            value = payload[key]
            self.assertTrue(value is None or (type(value) in (int, float) and math.isfinite(value)), key)
            if value is not None:
                self.assertGreaterEqual(value, 0, key)
        for key in ("dce_pct", "exclusion_rate_pct", "false_positive_rate_pct"):
            if payload[key] is not None:
                self.assertLessEqual(payload[key], 100, key)
        self.assertEqual(payload["chunks"], payload["unique_chunks"] + payload["duplicate_chunks"])
        self.assertEqual(
            payload["unique_chunks"],
            payload["delta_chunks"] + payload["self_compressed_chunks"] + payload["raw_chunks"],
        )
        self.assertLessEqual(payload["similar_chunks"], payload["unique_chunks"])
        self.assertLessEqual(payload["delta_chunks"], payload["similar_chunks"])
        self.assertEqual(payload["input_bytes"], payload["source_read_bytes"])
        self.assertLessEqual(payload["after_dedup_bytes"], payload["input_bytes"])
        self.assertLessEqual(payload["after_delta_bytes"], payload["after_dedup_bytes"])
        self.assertLessEqual(payload["estimated_final_payload_bytes"], payload["after_dedup_bytes"])
        final_bytes = payload["estimated_final_payload_bytes"]
        if final_bytes:
            self.assertAlmostEqual(payload["input_bytes"] / final_bytes, payload["estimated_drr"], places=9)
        else:
            self.assertIsNone(payload["estimated_drr"])
        after_dedup = payload["after_dedup_bytes"]
        if after_dedup:
            expected_dce = (1 - payload["after_delta_bytes"] / after_dedup) * 100
            self.assertAlmostEqual(expected_dce, payload["dce_pct"], places=9)
        else:
            self.assertIsNone(payload["dce_pct"])
        expected_checks = (
            payload["delta_chunks"] + payload["self_compressed_chunks"] if payload["verify"] else 0
        )
        self.assertEqual(expected_checks, payload["verification_checks"])
        self.assertEqual(0, payload["verification_failures"])

    def test_help_and_argument_exit_codes(self):
        sample = self.fixture("arguments.bin", b"x")
        self.assertIn("usage", self.invoke("--help").stdout.lower())
        invalid = (
            ((), 2),
            (("--input", sample, "--scheme", "unknown"), 2),
            (("--input", sample, "--scheme", "od-g"), 2),
            (("--input", sample, "--scheme", "od-g", "--features", "-1", "--super-features", "1"), 2),
            (("--input", sample, "--scheme", "od-g", "--features", "10", "--super-features", "3"), 2),
            (("--input", sample, "--scheme", "ss-g", *OD_ARGS), 2),
            (("--input", self.tmp / "missing.bin", "--scheme", "ss-g"), 1),
        )
        for args, status in invalid:
            with self.subTest(args=args):
                self.invoke(*args, expected=status)

    def test_json_schema_and_default_verification(self):
        path = self.fixture('schema-"quoted\\tab\t.bin', repeating(128 * 1024))
        payload = self.run_json(path, "ss-g-s", pipeline=None)
        self.assertEqual(str(path), payload["input"])
        self.assertEqual("parallel", payload["pipeline"])
        self.assertTrue(payload["verify"])
        self.assertGreater(payload["self_compressed_chunks"], 0)
        self.assertIsNone(payload["features"])
        self.assertIsNone(payload["super_features"])

    def test_no_verify(self):
        path = self.fixture("no-verify.bin", repeating(128 * 1024))
        payload = self.run_json(path, "ss-g-s", no_verify=True)
        self.assertFalse(payload["verify"])
        self.assertEqual(0, payload["verification_checks"])

    def test_empty_input(self):
        payload = self.run_json(self.fixture("empty.bin", b""))
        for key in COUNT_KEYS:
            self.assertEqual(0, payload[key], key)
        for key in RATE_KEYS:
            self.assertIsNone(payload[key], key)

    def test_size_boundaries(self):
        kib = 1024
        mib = 1024 * kib
        sizes = (1, 4 * kib - 1, 4 * kib, 4 * kib + 1, 8 * kib - 1, 8 * kib, 8 * kib + 1,
                 32 * kib - 1, 32 * kib, 32 * kib + 1, 4 * mib - 1, 4 * mib, 4 * mib + 1)
        for size in sizes:
            with self.subTest(size=size):
                path = self.fixture(f"zero-{size}.bin", bytes(size))
                self.assertEqual(size, self.run_json(path)["input_bytes"])

    def test_content_patterns(self):
        size = 256 * 1024
        patterns = {"zero": bytes(size), "repeating": repeating(size), "fixed-random": fixed_random(size)}
        for name, data in patterns.items():
            with self.subTest(pattern=name):
                path = self.fixture(f"{name}.bin", data)
                self.assertEqual(size, self.run_json(path, "ss-g-s", "parallel")["input_bytes"])

    def test_multiple_input_buffers(self):
        size = 8 * 1024 * 1024 + 257
        path = self.fixture("multi-buffer.bin", repeating(size))
        self.assertEqual(size, self.run_json(path, pipeline="parallel")["input_bytes"])

    def test_all_schemes_match_between_pipelines(self):
        block = bytearray(fixed_random(512 * 1024))
        data = bytearray(block * 4)
        for copy in range(1, 4):
            for offset in (4093, 65539, 262147, 400009):
                data[copy * len(block) + offset] ^= copy
        path = self.fixture("pipeline-equivalence.bin", data)
        for scheme in SCHEMES:
            with self.subTest(scheme=scheme):
                sequential = self.run_json(path, scheme, "sequential")
                parallel = self.run_json(path, scheme, "parallel")
                self.assertGreater(sequential["delta_chunks"], 0)
                logical_seq = {key: value for key, value in sequential.items() if key not in TIMING_KEYS | {"pipeline"}}
                logical_parallel = {key: value for key, value in parallel.items() if key not in TIMING_KEYS | {"pipeline"}}
                self.assertEqual(logical_seq, logical_parallel)
                if scheme.startswith("od-"):
                    self.assertEqual((12, 3), (sequential["features"], sequential["super_features"]))


def main():
    global BINARY, TIMEOUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="./speedsketch", help="SpeedSketch CLI path (default: ./speedsketch)")
    parser.add_argument("--timeout", type=float, default=30.0, help="timeout per CLI invocation in seconds")
    args, unittest_args = parser.parse_known_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    BINARY = Path(args.binary).expanduser().resolve()
    TIMEOUT = args.timeout
    unittest.main(argv=[sys.argv[0], *unittest_args], verbosity=2)


if __name__ == "__main__":
    main()
