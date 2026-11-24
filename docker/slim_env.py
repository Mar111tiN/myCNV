"""Slim a conda env for a runtime image: drop files nothing runs, strip ELF debug symbols (plan 3b, 2026-11-24).

    python -I slim_env.py <prefix> [--remove PKG ...] [--keep-pycache] [--dry-run]

--remove: conda packages to force-remove first (`micromamba remove -y -p <prefix> --force`), only those that are
installed (a changed solver result never fails the build; missing names are logged). Use it for packages nothing
loads at runtime: pytest and its helpers, tk (no tkinter use), cairo with its X11 and font packages (myCNV: only the
optional cairo build of harfbuzz needs them; matplotlib links raqm/harfbuzz/glib, which stay).

Removed: the package cache pkgs/ (hardlinked with the env; not needed at runtime), test suites of the packages and of
the stdlib (site-packages/**/tests, lib/python3.*/test), idlelib,
turtledemo, tkinter, ensurepip, C headers (include/), static libraries (*.a), docs/man/info, translations
(share/locale), *.js.map, and __pycache__ (unless --keep-pycache; the image then sets PYTHONDONTWRITEBYTECODE=1:
about 1 s more per Python start, measured on the workbench). Every ELF file gets `strip --strip-unneeded`, which
keeps the dynamic symbol table, so dlopen/dlsym and ctypes still work. conda-meta is kept (provenance, Trivy).
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys


def tree_size(path: str) -> int:
    """Bytes under path, every inode once (micromamba hardlinks the env files with its package cache pkgs/)."""
    total, seen = 0, set()
    for root, _, files in os.walk(path):
        for name in files:
            p = os.path.join(root, name)
            if os.path.islink(p):
                continue
            st = os.lstat(p)
            if (st.st_dev, st.st_ino) in seen:
                continue
            seen.add((st.st_dev, st.st_ino))
            total += st.st_size
    return total


def is_elf(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(4) == b"\x7fELF"
    except OSError:
        return False


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("prefix")
    ap.add_argument("--remove", nargs="+", default=[], metavar="PKG")
    ap.add_argument("--keep-pycache", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    p = os.path.abspath(a.prefix)
    if a.remove:
        meta = os.path.join(p, "conda-meta")
        installed = {"-".join(f[:-5].split("-")[:-2]) for f in os.listdir(meta) if f.endswith(".json")}
        present = [pkg for pkg in a.remove if pkg in installed]
        missing = [pkg for pkg in a.remove if pkg not in installed]
        if missing:
            print(f"slim_env: not installed, skipped: {' '.join(missing)}", file=sys.stderr)
        if present and not a.dry_run:
            subprocess.run(["micromamba", "remove", "-y", "-p", p, "--force", *present], check=True)
        print(f"slim_env: removed packages {' '.join(present) or '-'}", file=sys.stderr)
    stdlib = sorted(glob.glob(os.path.join(p, "lib", "python3.*")))
    if not stdlib:
        raise SystemExit(f"{p}: no lib/python3.* (not a conda env with Python?)")
    py = stdlib[-1]
    site = os.path.join(py, "site-packages")
    before = tree_size(p)

    # pkgs/ = micromamba's package cache, hardlinked with the env files. Without deleting it, every stripped or
    # deleted env file stays behind there as a second copy (CI 2026-11-24: mycnv image 873 MB -> 1.02 GB).
    dirs = [os.path.join(p, "pkgs"), os.path.join(p, "include"), os.path.join(p, "man")]
    dirs += [os.path.join(p, "share", d) for d in ("doc", "man", "info", "gtk-doc", "locale")]
    dirs += [os.path.join(py, d) for d in ("test", "idlelib", "turtledemo", "tkinter", "ensurepip")]
    for root, subdirs, _ in os.walk(site):
        for d in list(subdirs):
            if d == "tests" or (d == "__pycache__" and not a.keep_pycache):
                dirs.append(os.path.join(root, d))
                subdirs.remove(d)
    if not a.keep_pycache:
        for root, subdirs, _ in os.walk(p):
            if root.startswith(site) or "conda-meta" in root:
                continue
            for d in list(subdirs):
                if d == "__pycache__":
                    dirs.append(os.path.join(root, d))
                    subdirs.remove(d)
    files = []
    for root, _, names in os.walk(p):
        if "conda-meta" in root:
            continue
        for n in names:
            if n.endswith((".a", ".js.map")) or (n.endswith(".pyc") and not a.keep_pycache):
                files.append(os.path.join(root, n))
    removed = 0
    for d in sorted(set(dirs)):
        if os.path.isdir(d) and not os.path.islink(d):
            removed += tree_size(d)
            if not a.dry_run:
                shutil.rmtree(d)
    for f in files:
        if os.path.lexists(f):
            removed += 0 if os.path.islink(f) else os.path.getsize(f)
            if not a.dry_run:
                os.remove(f)

    stripped = saved = 0
    failed = []
    for root, _, names in os.walk(p):
        if "conda-meta" in root:
            continue
        for n in names:
            f = os.path.join(root, n)
            if os.path.islink(f) or not os.path.isfile(f) or not is_elf(f):
                continue
            size = os.path.getsize(f)
            if a.dry_run:
                continue
            r = subprocess.run(["strip", "--strip-unneeded", f], capture_output=True, text=True)
            if r.returncode:
                failed.append(f"{os.path.relpath(f, p)}: {r.stderr.strip()[:120]}")
                continue
            stripped += 1
            saved += size - os.path.getsize(f)
    after = tree_size(p)
    print(f"slim_env {p}: {before / 1e6:.0f} MB -> {after / 1e6:.0f} MB (removed files {removed / 1e6:.0f} MB, "
          f"strip saved {saved / 1e6:.0f} MB on {stripped} ELF files; {len(failed)} not strippable)", file=sys.stderr)
    for line in failed[:20]:
        print(f"  not stripped: {line}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
