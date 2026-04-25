import os
from pathlib import Path
import shutil
import subprocess


def test_gallery_launcher_rebuilds_index_with_structure_metrics(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    source_gallery = repo_root / "gallery"

    launcher_root = tmp_path / "launcher"
    launcher_root.mkdir()

    gallery_script = launcher_root / "gallery"
    shutil.copy2(source_gallery, gallery_script)
    gallery_script.chmod(0o755)

    data_dir = launcher_root / "data"
    data_dir.mkdir()

    molstar_dir = launcher_root / "node_modules" / "molstar" / "build" / "viewer"
    molstar_dir.mkdir(parents=True)
    (molstar_dir / "molstar.js").write_text("window.molstar = {};\n", encoding="utf-8")
    (molstar_dir / "molstar.css").write_text("body {}\n", encoding="utf-8")
    playwright_dir = launcher_root / "node_modules" / "playwright"
    playwright_dir.mkdir(parents=True)
    (playwright_dir / "package.json").write_text("{}", encoding="utf-8")

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    summary_path = run_dir / "run.summary.json"
    summary_path.write_text("{}", encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_uv = bin_dir / "uv"
    fake_node = bin_dir / "node"
    log_path = tmp_path / "uv.log"
    fake_uv.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "printf '%s\\n' \"$*\" >> \"$FAKE_UV_LOG\"\n",
        encoding="utf-8",
    )
    fake_uv.chmod(0o755)
    fake_node.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "exit 0\n",
        encoding="utf-8",
    )
    fake_node.chmod(0o755)

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["FAKE_UV_LOG"] = str(log_path)

    subprocess.run(
        [
            "bash",
            str(gallery_script),
            str(summary_path),
            "--force-index",
            "--skip-thumbnails",
        ],
        check=True,
        cwd=launcher_root,
        env=env,
        capture_output=True,
        text=True,
    )

    invocations = log_path.read_text(encoding="utf-8").splitlines()
    index_call = next(
        invocation
        for invocation in invocations
        if " -m volumizer.gallery_index " in f" {invocation} "
    )

    assert "--summary" in index_call
    assert "--db" in index_call
    assert "--progress-every 500" in index_call
    assert "--skip-structure-metrics" not in index_call


def test_gallery_launcher_skip_thumbnails_relinks_existing_cached_outputs(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    source_gallery = repo_root / "gallery"

    launcher_root = tmp_path / "launcher"
    launcher_root.mkdir()

    gallery_script = launcher_root / "gallery"
    shutil.copy2(source_gallery, gallery_script)
    gallery_script.chmod(0o755)

    data_dir = launcher_root / "data"
    data_dir.mkdir()

    molstar_dir = launcher_root / "node_modules" / "molstar" / "build" / "viewer"
    molstar_dir.mkdir(parents=True)
    (molstar_dir / "molstar.js").write_text("window.molstar = {};\n", encoding="utf-8")
    (molstar_dir / "molstar.css").write_text("body {}\n", encoding="utf-8")
    playwright_dir = launcher_root / "node_modules" / "playwright"
    playwright_dir.mkdir(parents=True)
    (playwright_dir / "package.json").write_text("{}", encoding="utf-8")

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    summary_path = run_dir / "run.summary.json"
    summary_path.write_text("{}", encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_uv = bin_dir / "uv"
    fake_node = bin_dir / "node"
    log_path = tmp_path / "uv.log"
    fake_uv.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "printf '%s\\n' \"$*\" >> \"$FAKE_UV_LOG\"\n",
        encoding="utf-8",
    )
    fake_uv.chmod(0o755)
    fake_node.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "exit 0\n",
        encoding="utf-8",
    )
    fake_node.chmod(0o755)

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["FAKE_UV_LOG"] = str(log_path)

    subprocess.run(
        [
            "bash",
            str(gallery_script),
            str(summary_path),
            "--force-index",
            "--skip-thumbnails",
        ],
        check=True,
        cwd=launcher_root,
        env=env,
        capture_output=True,
        text=True,
    )

    invocations = log_path.read_text(encoding="utf-8").splitlines()
    render_call = next(
        invocation
        for invocation in invocations
        if " -m volumizer.gallery_render " in f" {invocation} "
    )

    assert "--db" in render_call
    assert "--reuse-existing-only" in render_call
    assert "--progress-interval 0" in render_call


def test_gallery_launcher_include_failed_forwards_to_renderer(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    source_gallery = repo_root / "gallery"

    launcher_root = tmp_path / "launcher"
    launcher_root.mkdir()

    gallery_script = launcher_root / "gallery"
    shutil.copy2(source_gallery, gallery_script)
    gallery_script.chmod(0o755)

    data_dir = launcher_root / "data"
    data_dir.mkdir()

    molstar_dir = launcher_root / "node_modules" / "molstar" / "build" / "viewer"
    molstar_dir.mkdir(parents=True)
    (molstar_dir / "molstar.js").write_text("window.molstar = {};\n", encoding="utf-8")
    (molstar_dir / "molstar.css").write_text("body {}\n", encoding="utf-8")
    playwright_dir = launcher_root / "node_modules" / "playwright"
    playwright_dir.mkdir(parents=True)
    (playwright_dir / "package.json").write_text("{}", encoding="utf-8")

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    summary_path = run_dir / "run.summary.json"
    summary_path.write_text("{}", encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_uv = bin_dir / "uv"
    fake_node = bin_dir / "node"
    log_path = tmp_path / "uv.log"
    fake_uv.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "printf '%s\\n' \"$*\" >> \"$FAKE_UV_LOG\"\n",
        encoding="utf-8",
    )
    fake_uv.chmod(0o755)
    fake_node.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "exit 0\n",
        encoding="utf-8",
    )
    fake_node.chmod(0o755)

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["FAKE_UV_LOG"] = str(log_path)

    subprocess.run(
        [
            "bash",
            str(gallery_script),
            str(summary_path),
            "--force-index",
            "--include-failed",
        ],
        check=True,
        cwd=launcher_root,
        env=env,
        capture_output=True,
        text=True,
    )

    invocations = log_path.read_text(encoding="utf-8").splitlines()
    render_call = next(
        invocation
        for invocation in invocations
        if " -m volumizer.gallery_render " in f" {invocation} "
    )

    assert "--include-failed" in render_call
