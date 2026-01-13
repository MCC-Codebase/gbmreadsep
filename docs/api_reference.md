# Programmatic API reference

This page documents public Python entrypoints useful for scripting.

---

## `gbmreadsep.toy_data.make_toy_data(outdir)`

**Purpose**: Generate a tiny reference, BAM, and VCF for demos/tests.

**Parameters**
- `outdir` (str | Path): Output directory.

**Returns**
- `dict` with paths to the generated files.

**Errors**
- Raises `OSError` if output paths cannot be created.

**Example**
```python
from gbmreadsep.toy_data import make_toy_data

summary = make_toy_data(outdir="toy")
print(summary["anchors_vcf"])
```

**Source Anchor**: src/gbmreadsep/toy_data.py:L45-L149#b15126af

---

## `gbmreadsep.nanopore.align.build_alignment_commands(...)`

**Purpose**: Return minimap2/samtools commands without executing them.

**Parameters**
- `fastq`: list of FASTQ paths
- `ref_fa`: reference FASTA
- `out_bam`: output BAM path
- `threads`, `minimap2_preset`, `sort_mem`

**Returns**
- `dict` with `minimap2`, `samtools_sort`, and `samtools_index` command lists.

**Errors**
- Raises `ValueError` on unsupported presets.

**Example**
```python
from gbmreadsep.nanopore.align import build_alignment_commands

cmds = build_alignment_commands(
    fastq=["reads.fq.gz"],
    ref_fa="ref.fa",
    out_bam="tumor.bam",
)
print(cmds["minimap2"])
```

**Source Anchor**: src/gbmreadsep/nanopore/align.py:L35-L95#614a8059

---

## `gbmreadsep.nanopore.somatic.build_somatic_commands(...)`

**Purpose**: Return somatic caller commands without executing them.

**Parameters**
- `tumor_bam`, `normal_bam` (optional)
- `ref_fa`, `outdir`
- `engine` (`docker` or `native`)

**Returns**
- `dict` including either `docker_cmd` or `cmd` plus expected output paths.

**Errors**
- Raises `FileNotFoundError` for missing input paths.

**Example**
```python
from gbmreadsep.nanopore.somatic import build_somatic_commands

cmds = build_somatic_commands(
    tumor_bam="tumor.bam",
    ref_fa="ref.fa",
    outdir="somatic_vcf",
    engine="docker",
)
print(cmds["docker_cmd"])
```

**Source Anchor**: src/gbmreadsep/nanopore/somatic.py:L96-L265#ba7a0962

---

## `gbmreadsep.anchor_index.build_index(...)` and `verify_index(...)`

**Purpose**: Generate and validate source anchors for public symbols.

**Parameters**
- `build_index(paths, out_path, include_private=False, include_methods=True)`
- `verify_index(index_path)`

**Returns**
- `build_index`: mapping of symbol → metadata (and writes JSON)
- `verify_index`: report dict; raises `SystemExit(1)` on mismatches

**Example**
```python
from gbmreadsep.anchor_index import build_index, verify_index

build_index(["src/gbmreadsep"], "docs/_generated/anchors_index.json")
verify_index("docs/_generated/anchors_index.json")
```

**Source Anchors**:
- build_index: src/gbmreadsep/anchor_index.py:L52-L119#29fbc4d7
- verify_index: src/gbmreadsep/anchor_index.py:L122-L150#450e2586
