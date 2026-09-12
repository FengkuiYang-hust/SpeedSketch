# SpeedSketch

SpeedSketch 是论文 *SpeedSketch: An Ultra-Fast Sketch Generation and Delta Encoding Framework for Delta Compression* 的开源研究实现。当前命令行程序用于在单个输入文件上比较 sketch/candidate search 与 Delta Encoding 组合。

当前 FastCDC 的 min/avg/max 固定为 `4/8/32 KiB`。`4 KiB` 最小块是工程折中，有意偏离论文使用的 `2/8/32 KiB`；因此默认结果不声称是论文原始 FastCDC 参数下的复现。

## 构建与测试

需要 GNU Make、支持 C++14/C99 的编译器和 zstd 开发库。

```sh
make -j
make check
```

产物为 `./speedsketch`。CLI 黑盒测试只依赖 Python 标准库，也可直接指定 binary 和每次调用的超时：

```sh
python3 tests/test_cli.py --binary ./speedsketch --timeout 30
```

`make asan` 和 `make tsan` 分别运行 Address/UndefinedBehavior Sanitizer 与 Thread Sanitizer 检查；TSan target 使用 Linux `setarch` 关闭 ASLR，以避开本机 libtsan 的地址映射冲突。

## CLI

```text
speedsketch --input PATH --scheme SCHEME
            [--pipeline sequential|parallel]
            [--features N --super-features N]
            [--json] [--no-verify]
```

- `--input` 与 `--scheme` 必填；`--pipeline` 默认为 `parallel`。
- `--json` 令 stdout 只输出一个 JSON object；诊断信息写 stderr。
- 默认会对每个最终选中的 delta 或 self-compressed representation 执行 decode/decompress，并与原 chunk 逐字节比较。`--no-verify` 会关闭该门禁，只应用于明确标注的探索性计时。
- `od-g`、`od-x` 必须显式提供正整数 `--features` 与 `--super-features`，且前者不小于并可整除后者；其他 scheme 拒绝这两个参数。
- 正常完成返回 `0`，参数/配置错误返回 `2`，I/O、运行或验证失败返回 `1`。
- 输入必须在一次运行中保持不变；程序会校验最终读取量，并拒绝无法补齐的读取或提前 EOF，但不会为同尺寸并发覆写创建快照。

Scheme 与论文组合一致：

| Scheme | 相似候选 | Delta encoder |
|---|---|---|
| `ss-g` | SpeedSketch | native Gdelta |
| `ss-g-s` | SpeedSketch | sketch-accelerated Gdelta |
| `ss-x` | SpeedSketch | Xdelta |
| `od-g` | Odess | native Gdelta |
| `od-x` | Odess | Xdelta |

示例：

```sh
./speedsketch --input data.bin --scheme ss-g-s --pipeline parallel --json
./speedsketch --input data.bin --scheme od-g --features 12 --super-features 3 --json
```

测试采用 `features=12`、`super-features=3` 作为一个合法、固定的 OD 覆盖配置；这不是论文声明的唯一参数，也不代表其他参数已经完成质量或性能验证。

## JSON metrics

所有 schema key 始终存在；不适用或分母为零的值为 JSON `null`。

- 身份与配置：`schema_version`、`input`、`scheme`、`pipeline`、`verify`、`features`、`super_features`。
- 输入与分类：`input_bytes`、`chunks`、`unique_chunks`、`duplicate_chunks`、`similar_chunks`、`delta_chunks`、`self_compressed_chunks`、`raw_chunks`。`similar_chunks` 是找到 sketch/SF base 的 unique chunk 数；`delta_chunks` 只统计 delta payload 严格优于 raw 和 self-compressed 后被最终选中的数量。
- 字节量：`after_dedup_bytes` 是 unique chunks 的原始字节和；`after_delta_bytes` 对有 base 的 chunk 累计 `min(raw chunk, raw delta)`，无 base 时累计 raw，因此不混入 ZSTD self-compression。最终估算对 G schemes 在 raw、`ZSTD(raw)`、`ZSTD(raw Gdelta)` 中择小，对 X schemes 在 raw、`ZSTD(raw)`、raw Xdelta 中择小；duplicates 计零。`estimated_final_payload_bytes` 不含持久化索引、recipe、container padding、journal 或其他 Archive metadata。
- 读取量：`source_read_bytes` 是 stream reader 从输入顺序读取的字节，正常完成时等于 `input_bytes`；`base_read_bytes` 累计为 XXH64 collision/duplicate 逐字节确认和 similar-base encode 所做的 `pread`。
- 缩减效果：`estimated_drr = input_bytes / estimated_final_payload_bytes`；`dce_pct = (1 - after_delta_bytes / after_dedup_bytes) * 100`。
- lookup 过滤：`exclusion_rate_pct = excluded_lookups / lookup_judgments * 100`；`false_positive_rate_pct = false_positive_lookups / lookup_judgments * 100`，与论文式 (15) 一样以全部 judgments 为分母。
- 时间与吞吐：`wall_seconds`、`chunking_seconds`、`dedup_seconds`、`sketch_seconds`、`encoding_seconds`、`verification_seconds`；wall 与前四个 pipeline stage 另有对应的 `*_mib_s`。这些字段会受 cache、调度和 verify 设置影响。
- 正确性：`verification_checks` 与 `verification_failures`；默认验证时 checks 等于选中的 delta 与 self-compressed chunks 之和，正常运行要求 failures 为零。

`sequential` 与 `parallel` 使用相同算法和顺序语义；除 pipeline 名称、时间和吞吐字段外，逻辑 metrics 应完全一致。

## 可行性输入与结果边界

下面的脚本生成固定的 64 MiB 重复数据并施加稀疏、确定性修改，同时打印 SHA256；fixture 不写入仓库，且脚本拒绝覆盖已有路径：

```sh
python3 tests/generate_feasibility.py /tmp/speedsketch-feasibility-64MiB.bin
# 94da0c0d90d1728fd924feac13100afada1a79a18baa4a6b0243ac6820a80f01

./speedsketch --input /tmp/speedsketch-feasibility-64MiB.bin \
  --scheme ss-g-s --pipeline parallel --json
```

这只构成 `0 warmup + 1 measure` 的可行性观测：保留单次原始输出，不计算平均值、置信区间，也不外推到论文完整数据集。正式性能结论需要固定环境和输入、基线、至少 `1 warmup + 5 measures`、原始结果及逐字节验证。

当前程序是算法 simulator，不是可持久化的 Archive/restore 系统：它不实现 durable container、recipe/index 落盘、崩溃一致性、GC 或独立 restore 路径。因此 estimated payload、DRR、读字节和阶段吞吐不能解释为最终物理 Archive、恢复 I/O 或端到端 durable performance。
