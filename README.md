# SpeedSketch

SpeedSketch is the open-source research implementation of *SpeedSketch: An Ultra-Fast Sketch Generation and Delta Encoding Framework for Delta Compression*. The current command-line program compares combinations of sketch or candidate search methods and delta encoders on a single input file.

FastCDC currently uses fixed minimum, average, and maximum chunk sizes of `4/8/32 KiB`. The `4 KiB` minimum is an engineering tradeoff and intentionally differs from the paper's `2/8/32 KiB` configuration. Results produced with the default settings should therefore not be presented as a reproduction of the paper's original FastCDC configuration.

## Build and Test

The build requires GNU Make, compilers with C++14 and C99 support, and the zstd development library.

```sh
make -j
make check
```

The resulting executable is `./speedsketch`. The CLI black-box tests depend only on the Python standard library. You can also specify the executable and per-command timeout directly:

```sh
python3 tests/test_cli.py --binary ./speedsketch --timeout 30
```

`make asan` runs the AddressSanitizer and UndefinedBehaviorSanitizer checks, while `make tsan` runs the ThreadSanitizer checks. The TSan target uses Linux `setarch` to disable ASLR and avoid a libtsan address-mapping conflict on this host.

## CLI

```text
speedsketch --input PATH --scheme SCHEME
            [--pipeline sequential|parallel]
            [--features N --super-features N]
            [--json] [--no-verify]
```

- `--input` and `--scheme` are required. `--pipeline` defaults to `parallel`.
- `--json` makes stdout contain exactly one JSON object. Diagnostics are written to stderr.
- By default, the program decodes or decompresses every selected delta or self-compressed representation and compares it byte for byte with the original chunk. `--no-verify` disables this gate and should be used only for explicitly labeled exploratory timing runs.
- `od-g` and `od-x` require positive integer values for `--features` and `--super-features`. The feature count must be at least the super-feature count and must be evenly divisible by it. Other schemes reject these options.
- A successful run returns `0`, an argument or configuration error returns `2`, and an I/O, runtime, or verification failure returns `1`.
- The input must remain unchanged throughout a run. The program validates the final byte count and rejects reads that cannot be completed or encounter an early EOF, but it does not snapshot concurrent same-size overwrites.

The schemes match the combinations evaluated in the paper:

| Scheme | Similarity candidate search | Delta encoder |
|---|---|---|
| `ss-g` | SpeedSketch | Native Gdelta |
| `ss-g-s` | SpeedSketch | Sketch-accelerated Gdelta |
| `ss-x` | SpeedSketch | Xdelta |
| `od-g` | Odess | Native Gdelta |
| `od-x` | Odess | Xdelta |

Examples:

```sh
./speedsketch --input data.bin --scheme ss-g-s --pipeline parallel --json
./speedsketch --input data.bin --scheme od-g --features 12 --super-features 3 --json
```

The tests use `features=12` and `super-features=3` as one valid, fixed Odess configuration. This is not claimed to be the paper's only configuration, and other values have not necessarily passed the same correctness or performance validation.

## JSON Metrics

Every schema key is always present. Values that do not apply or have a zero denominator are JSON `null`.

- Identity and configuration: `schema_version`, `input`, `scheme`, `pipeline`, `verify`, `features`, and `super_features`.
- Input and classification: `input_bytes`, `chunks`, `unique_chunks`, `duplicate_chunks`, `similar_chunks`, `delta_chunks`, `self_compressed_chunks`, and `raw_chunks`. `similar_chunks` counts unique chunks for which a sketch or super-feature base was found. `delta_chunks` counts only delta payloads that were ultimately selected because they were strictly smaller than both the raw and self-compressed representations.
- Byte counts: `after_dedup_bytes` is the sum of the original sizes of all unique chunks. For a chunk with a base, `after_delta_bytes` adds `min(raw chunk, raw delta)`; without a base, it adds the raw chunk size. This metric therefore excludes ZSTD self-compression. The final estimate for Gdelta schemes selects the smallest of raw data, `ZSTD(raw)`, and `ZSTD(Gdelta(raw))`; Xdelta schemes select the smallest of raw data, `ZSTD(raw)`, and raw Xdelta output. Duplicates contribute zero bytes. `estimated_final_payload_bytes` excludes persistent indexes, recipes, container padding, journals, and all other archive metadata.
- Read volume: `source_read_bytes` is the number of bytes read sequentially from the input and equals `input_bytes` after a successful run. `base_read_bytes` counts bytes read with `pread` to confirm XXH64 collision or duplicate candidates and to encode against similar bases.
- Data reduction: `estimated_drr = input_bytes / estimated_final_payload_bytes`; `dce_pct = (1 - after_delta_bytes / after_dedup_bytes) * 100`.
- Lookup filtering: `exclusion_rate_pct = excluded_lookups / lookup_judgments * 100`; `false_positive_rate_pct = false_positive_lookups / lookup_judgments * 100`. Like equation (15) in the paper, the latter uses all judgments as its denominator.
- Timing and throughput: `wall_seconds`, `chunking_seconds`, `dedup_seconds`, `sketch_seconds`, `encoding_seconds`, and `verification_seconds`. The wall time and first four pipeline stages also have corresponding `*_mib_s` fields. These values are affected by cache state, scheduling, and verification settings.
- Correctness: `verification_checks` and `verification_failures`. With default verification, the number of checks equals the total number of selected delta and self-compressed chunks, and a successful run requires zero failures.

The `sequential` and `parallel` pipelines use the same algorithms and ordering semantics. Except for the pipeline name, timing, and throughput fields, their logical metrics should be identical.

## Feasibility Fixture and Result Scope

The following script generates a fixed 64 MiB data set with repeated regions and sparse deterministic modifications, then prints its SHA256 digest. The fixture is not stored in the repository, and the script refuses to overwrite an existing path.

```sh
python3 tests/generate_feasibility.py /tmp/speedsketch-feasibility-64MiB.bin
# 94da0c0d90d1728fd924feac13100afada1a79a18baa4a6b0243ac6820a80f01

./speedsketch --input /tmp/speedsketch-feasibility-64MiB.bin \
  --scheme ss-g-s --pipeline parallel --json
```

This is only a `0 warmup + 1 measure` feasibility observation. Preserve the raw output from the single run; do not calculate an average or confidence interval, and do not extrapolate it to the paper's full data set. A formal performance result requires a fixed environment and input, a baseline, at least `1 warmup + 5 measures`, raw results, and byte-for-byte verification.

The current program is an algorithm simulator, not a persistent archive and restore system. It does not implement durable containers, persistent recipes or indexes, crash consistency, garbage collection, or an independent restore path. The estimated payload, DRR, read volume, and stage throughput must not be interpreted as final physical archive size, restore I/O, or end-to-end durable performance.
