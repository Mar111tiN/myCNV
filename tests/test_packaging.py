"""Image and packaging contract: no patient data in the build context, lean imports, pinned stages."""

import os
import re
import subprocess
import sys

import pytest

from conftest import HERE

ROOT = os.path.dirname(HERE)

# The Dockerfile and .dockerignore are checked in the source checkout (CI pytest job, workbench). The image's test
# stage copies only src/ and tests/ (.dockerignore excludes both files), so these checks skip there.
needs_checkout = pytest.mark.skipif(
    not (os.path.isfile(os.path.join(ROOT, "Dockerfile")) and os.path.isfile(os.path.join(ROOT, ".dockerignore"))),
    reason="Dockerfile/.dockerignore are not in this tree (Docker test stage)",
)


def _lines(name):
    with open(os.path.join(ROOT, name)) as fh:
        return [line.strip() for line in fh if line.strip() and not line.lstrip().startswith("#")]


@needs_checkout
def test_dockerignore_excludes_restricted_fixtures():
    assert "tests/data/restricted" in _lines(".dockerignore")


@needs_checkout
def test_dockerfile_copies_explicit_paths_and_runtime_is_last():
    lines = _lines("Dockerfile")
    copies = [line for line in lines if line.startswith("COPY") and "--from=" not in line]
    assert copies, "no COPY lines"
    for line in copies:
        assert not re.search(r"\s\.\s", line + " "), f"whole build context copied: {line}"
    stages = [line for line in lines if line.startswith("FROM")]
    assert stages[-1].endswith("AS runtime")
    for line in stages:
        image = line.split()[1]
        assert image == "build" or re.search(r":[0-9][\w.\-]*$", image), f"unpinned base image: {line}"
    assert any("MPLCONFIGDIR" in line for line in lines)


def _modules_after(code: str) -> set:
    out = subprocess.run([sys.executable, "-c", code + "\nimport sys\nprint(' '.join(sorted(sys.modules)))"],
                         check=True, capture_output=True, text=True).stdout
    return set(out.split())


def test_cli_without_plots_does_not_import_matplotlib():
    mods = _modules_after("from mycnv import cli, pipeline\ncli.main(['config'])")
    assert "matplotlib" not in mods


def test_plot_modules_do_not_need_seaborn():
    mods = _modules_after("import matplotlib\nmatplotlib.use('Agg')\nimport mycnv.plot, mycnv.plot_gistic")
    assert "seaborn" not in mods and "scipy" not in mods
