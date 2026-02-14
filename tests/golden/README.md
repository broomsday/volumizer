# Golden Outputs

This directory stores baseline annotation outputs used for regression/parity checks.

Current fixtures:
- `4jpn.annotation.json`: expected annotation dataframe at `2.0 A` resolution.

Regeneration workflow:
1. Run `bash src/compile_c_libs.sh`.
2. Run `uv run --python 3.11 python scripts/benchmark.py --group medium --resolution 2.0 --output-json AGENTS/benchmark.latest.json`.
3. Re-run `volumizer.volumize_pdb()` for the target fixture at the same resolution.
4. Save the resulting dataframe JSON into this directory.
5. Run `uv run --python 3.11 pytest tests/test_parity_scaffold.py` and verify tolerances still pass.
