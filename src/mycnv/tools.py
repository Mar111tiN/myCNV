"""Shared helpers: logging, mawk tool lookup, shell pipelines, table I/O and chromosome ordering."""

from __future__ import annotations

import logging
import os
import re
import shlex
import subprocess
from importlib import resources
from io import StringIO
from pathlib import Path

import pandas as pd

log = logging.getLogger("mycnv")

PRIMARY_CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX"]


def setup_logging(verbose: bool = False) -> None:
    """Log to stderr; data goes to files."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


def mawk_tool(name: str) -> Path:
    """Path of a packaged mawk tool (src/mycnv/mawk/<name>.mawk)."""
    path = resources.files("mycnv").joinpath("mawk", f"{name}.mawk")
    if not path.is_file():
        raise FileNotFoundError(f"mawk tool {name} not found in the mycnv package")
    return Path(str(path))


def mawk(name: str) -> str:
    """Shell snippet that runs a packaged mawk tool with bash (the tools use bash-only argument parsing)."""
    return f"bash {shlex.quote(str(mawk_tool(name)))}"


def _bash(cmd: str, **kwargs) -> subprocess.CompletedProcess:
    log.debug("$ %s", cmd)
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, **kwargs)


def run_cmd(cmd: str) -> None:
    """Run a shell pipeline with bash -o pipefail; raises CalledProcessError on failure."""
    _bash(cmd)


def cmd2df(cmd: str, **read_kwargs) -> pd.DataFrame:
    """Run a shell pipeline and read its tab-separated stdout (with header) into a DataFrame."""
    out = _bash(cmd, stdout=subprocess.PIPE).stdout.decode("utf-8")
    return pd.read_csv(StringIO(out), sep="\t", **read_kwargs)


def cat_cmd(path: str | os.PathLike) -> str:
    """Shell snippet that streams a plain or gzipped file to stdout."""
    path = str(path)
    return f"gunzip -c {shlex.quote(path)}" if path.endswith(".gz") else f"cat {shlex.quote(path)}"


def read_tsv(path: str | os.PathLike, **kwargs) -> pd.DataFrame:
    """Read a tab-separated table (plain or gzipped; pandas infers the compression)."""
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(df: pd.DataFrame, path: str | os.PathLike, **kwargs) -> None:
    """Write a tab-separated table without index; gzip if the name ends with .gz."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, **kwargs)


def chrom_key(chrom: str) -> tuple:
    """Natural sort key: chr1 < chr2 < ... < chr22 < chrX < chrY < chrM < others."""
    c = re.sub(r"^chr", "", str(chrom))
    if c.isdigit():
        return (0, int(c), "")
    order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    return (order.get(c, 4), 0, c)


def natural_chroms(chroms) -> list[str]:
    """Unique chromosome names in natural order."""
    return sorted({str(c) for c in chroms}, key=chrom_key)


def set_chrom_categorical(df: pd.DataFrame, chroms=None) -> pd.DataFrame:
    """Make Chr an ordered categorical (natural order) so sorting and groupby follow chr1..chr22, chrX.

    Uses column assignment, not df.loc[:, "Chr"] = ..., which keeps object dtype on pandas >= 2.
    """
    chroms = natural_chroms(df["Chr"].astype(str).unique()) if chroms is None else list(chroms)
    df = df.copy()
    df["Chr"] = pd.Categorical(df["Chr"].astype(str), categories=chroms, ordered=True)
    return df


def chrom_str(df: pd.DataFrame) -> pd.DataFrame:
    """Turn a categorical Chr back into plain strings (for writing and merging)."""
    if isinstance(df["Chr"].dtype, pd.CategoricalDtype):
        df = df.copy()
        df["Chr"] = df["Chr"].astype(str)
    return df
