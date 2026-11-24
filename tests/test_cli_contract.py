"""CLI contract (plan C1): the subcommands and flags that the kumo-variant-calling-nf modules call must not drift.

tests/data/cli_contract.json is the snapshot (made by tests/make_cli_contract.py from mycnv.cli.build_parser()).
"""

import argparse
import copy
import json
import shlex
import shutil
import types

import pytest

import mycnv
from make_cli_contract import (CONTRACT, NF_MODULES, SRC_PKG, TOP_LEVEL_MODULES, VERSION_MASK, affected_modules,
                               build_contract, contract_diff, render, source_note, stale_install)
from mycnv.cli import _check, build_parser


def mismatch_message(expected: dict, actual: dict, module=mycnv) -> str:
    lines = [f"the mycnv CLI no longer matches its contract {CONTRACT}:",
             *(f"  {c}" for c in contract_diff(expected, actual) or ["(see the rendered JSON)"]),
             *source_note(module)]
    if not stale_install(module):
        lines += [
            "If the change is deliberate: rerun `python tests/make_cli_contract.py` and commit the JSON with the change,",
            "and review the kumo-variant-calling-nf modules that call the changed flags (design section 3 process table):",
            *(f"  {m}" for m in affected_modules(expected, actual)),
        ]
    return "\n".join(lines)


def test_cli_contract_unchanged():
    expected = json.loads(CONTRACT.read_text())
    actual = build_contract()
    # as canonical JSON text, not `==`: in Python 1 == 1.0 == True, in the snapshot they differ
    assert render(actual) == render(expected), mismatch_message(expected, actual)


def test_contract_file_is_canonical():
    """The JSON is the script's own output (one line per option), so a rerun gives a minimal diff."""
    text = CONTRACT.read_text()
    assert text == render(json.loads(text)), "rerun `python tests/make_cli_contract.py` instead of editing the JSON"


def test_contract_covers_nested_commands_and_the_new_options():
    """`tracks gc|genmap` (nested) and the sex (M4/M5) and tumor-only (M1) options are part of the snapshot."""
    c = json.loads(CONTRACT.read_text())
    cmds = c["commands"]
    assert set(cmds) == {"ascat-input", "combine", "config", "pon-combine", "pon-coverage", "rolling", "sample", "sex",
                         "tn2cnv", "tracks"}
    assert set(cmds["tracks"]["commands"]) == {"gc", "genmap"} and cmds["tracks"]["command_required"] is True
    assert cmds["tracks"]["commands"]["genmap"]["options"]["--modes"]["default"] == "30_0,50_0,75_1,100_2"
    assert {"--cov", "--output"} <= set(cmds["sex"]["options"])
    assert "--out-sex" in cmds["combine"]["options"]
    for cmd in ("sample", "tn2cnv"):
        o = cmds[cmd]["options"]
        assert o["--tumor-only"]["action"] == "store_true" and o["--snp-sites"]["required"] is False
        assert o["--tumor-vaf"]["nargs"] == 2 and o["--tumor-vaf"]["type"] == "float"
    for cmd in ("combine", "rolling", "sex", "ascat-input"):
        assert "--tumor-only" in cmds[cmd]["options"] and "--snp-sites" not in cmds[cmd]["options"]
    assert cmds["pon-coverage"]["mutually_exclusive"] == [{"options": ["--matrix", "--ponmatrix"], "required": True}]


def iter_parsers(parser: argparse.ArgumentParser):
    """`parser` and all its subcommand parsers, nested ones included."""
    todo = [parser]
    while todo:
        p = todo.pop()
        yield p
        todo += [s for a in p._actions if isinstance(a, argparse._SubParsersAction) for s in a.choices.values()]


def exact_flags(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """`parser` and all its subcommand parsers without prefix matching (argparse's allow_abbrev).

    A module call must use the exact flags, else a stale flag keeps parsing after a rename to a longer one until a
    later option with the same prefix makes it ambiguous. mycnv itself turns prefix matching off (next test); this
    keeps the module-call check exact even if that changes.
    """
    for p in iter_parsers(parser):
        p.allow_abbrev = False
    return parser


def test_no_prefix_matching_anywhere():
    """mycnv's own rule (CHANGELOG: no argparse prefix matching), on every (nested) subcommand parser."""
    assert [p.prog for p in iter_parsers(build_parser()) if p.allow_abbrev] == []


# The module command lines, with placeholders for the Nextflow variables. The snapshot catches every change; these show
# directly that the module calls still parse and pass the CLI's cross-option checks (`_check`: e.g. --tumor-only needs
# --snp-sites). Sources: design section 3 (process table, 2026-11-23) for pon-coverage, pon-combine and the paired
# sample. The tumor-only sample line is NOT a section 3 line: it is the form that works after M1 (design 6.3 M1 plus
# `_check`): the paired line with --tumor-only --snp-sites added and --normal-name dropped (no --no-mock-germline
# either). Section 3 only notes "one-sample pileup + `--snp-sites $sites`", and that note read literally on the paired
# line does not parse (SECTION3_TUMOR_ONLY_READINGS below). Write the tumor-only module from this line until section 3
# spells it out.
MODULE_CALLS = {
    "MYCNV_PONCOVERAGE": "pon-coverage --matrix chr21.pon.gz --chrom chr21 --bed kit.padded.bed "
                         "--gc chr21.gc100-10.gz -o out -t 1",
    "MYCNV_PONCOMBINE": "pon-combine --ponmatrix pon -o out --chrom chr1 --chrom chr21 --chrom chrX",
    "MYCNV_SAMPLE paired": "sample --pileup P1.pileup.gz --ponmatrix ponmatrix --sample P1_T_vs_P1_N "
                           "--tumor-name P1_T --normal-name P1_N -o . --bed kit.padded.bed --gc gc --genmap genmap "
                           "--plots",
    "MYCNV_SAMPLE tumor-only": "sample --tumor-only --snp-sites kit.snp_sites.vcf.gz --pileup T1.pileup.gz "
                               "--ponmatrix ponmatrix --sample T1 --tumor-name T1 -o . --bed kit.padded.bed --gc gc "
                               "--genmap genmap --plots",
}


def parse_module_call(argv: list[str]) -> argparse.Namespace:
    parser = exact_flags(build_parser())
    args = parser.parse_args(argv)
    _check(parser, args)  # parser.error -> SystemExit, as in `mycnv` itself
    return args


@pytest.mark.parametrize("module", MODULE_CALLS)
def test_module_calls_parse(module):
    argv = shlex.split(MODULE_CALLS[module])
    assert parse_module_call(argv).cmd == argv[0]
    assert module.split()[0] == NF_MODULES[argv[0]].split()[0]  # the snapshot's module map names this module


def test_module_call_cross_checks_run():
    """The module-call check also applies the cross-option rules: the tumor-only call without its sites fails."""
    argv = shlex.split(MODULE_CALLS["MYCNV_SAMPLE tumor-only"].replace("--snp-sites kit.snp_sites.vcf.gz ", ""))
    with pytest.raises(SystemExit) as exit_:
        parse_module_call(argv)
    assert exit_.value.code == 2


def _flags(line: str) -> set[str]:
    return {t for t in shlex.split(line) if t.startswith("-")}


def test_tumor_only_call_is_the_paired_call_in_m1_form():
    """The tumor-only line is the paired section 3 line with --tumor-only --snp-sites added and --normal-name dropped.

    So the two module forms cannot drift apart unnoticed (a flag added to one only), and the comment on MODULE_CALLS
    says what the line is.
    """
    paired, tumor_only = MODULE_CALLS["MYCNV_SAMPLE paired"], MODULE_CALLS["MYCNV_SAMPLE tumor-only"]
    assert _flags(tumor_only) == (_flags(paired) - {"--normal-name"}) | {"--tumor-only", "--snp-sites"}
    assert "--no-mock-germline" not in _flags(tumor_only)


# Design section 3's tumor-only note ("one-sample pileup + `--snp-sites $sites`") read literally on the paired line:
# neither reading parses, so a module written from section 3 alone fails. The tumor-only line in MODULE_CALLS is the
# working form instead (test above).
SECTION3_TUMOR_ONLY_READINGS = {
    "paired + --snp-sites": ("--snp-sites kit.snp_sites.vcf.gz", "--snp-sites goes with --tumor-only"),
    "paired + --tumor-only --snp-sites": ("--tumor-only --snp-sites kit.snp_sites.vcf.gz",
                                          "--normal-name does not go with --tumor-only"),
}


@pytest.mark.parametrize("reading", SECTION3_TUMOR_ONLY_READINGS)
def test_section3_tumor_only_literal_readings_fail(reading, capsys):
    extra, error = SECTION3_TUMOR_ONLY_READINGS[reading]
    with pytest.raises(SystemExit) as exit_:
        parse_module_call(shlex.split(f"{MODULE_CALLS['MYCNV_SAMPLE paired']} {extra}"))
    assert exit_.value.code == 2
    assert error in capsys.readouterr().err


def test_version_output_matches_contract(capsys):
    """versions.yml parses `mycnv --version`: the printed line is the snapshot's template with the version."""
    template = json.loads(CONTRACT.read_text())["options"]["--version"]["version"]
    assert template == f"mycnv {VERSION_MASK}"
    with pytest.raises(SystemExit) as exit_:
        build_parser().parse_args(["--version"])
    assert exit_.value.code == 0
    assert capsys.readouterr().out == template.replace(VERSION_MASK, mycnv.__version__) + "\n"


def test_changes_name_the_calling_module():
    """A changed flag names the module of its subcommand; a nested one (`tracks genmap`) its top-level command."""
    old = build_contract()
    new = copy.deepcopy(old)
    sample = new["commands"]["sample"]
    del sample["options"]["--snp-sites"]
    sample["option_strings"].remove("--snp-sites")
    new["commands"]["tracks"]["commands"]["genmap"]["options"]["--modes"]["default"] = "30_0"
    assert contract_diff(old, new) == [
        "mycnv sample option_strings: removed --snp-sites",
        "mycnv sample --snp-sites: removed",
        'mycnv tracks genmap --modes default: "30_0,50_0,75_1,100_2" -> "30_0"',
    ]
    assert affected_modules(old, new) == [NF_MODULES["sample"],
                                          "`tracks`: no module calls it yet (genome-bundle track builds, plan D22)"]
    new = copy.deepcopy(old)
    new["commands"]["combine"]["options"]["--min-pon-match"]["default"] = 0.8
    assert affected_modules(old, new) == [
        "`combine`: no module calls it yet (step-by-step path; MYCNV_SAMPLE runs `sample`)"]


# ---------------------------------------------------------------- the generator itself, on a toy parser
def toy_parser(count_default=0, version="toy 1.2", tag_flag="--tag"):
    p = argparse.ArgumentParser(prog="toy")
    p.add_argument("--version", action="version", version=version)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("run")
    s.add_argument("infile")
    s.add_argument("-n", "--count", action="count", default=count_default)
    s.add_argument("--mode", action="store_const", const="fast")
    s.add_argument(tag_flag, dest="tag", action="extend", nargs="+", default=["a"])
    s.add_argument("--level", type=float, choices=[0.5, 1.0], default=0.5)
    s.add_argument("--clock", default=object())  # a computed default: not part of the contract
    g = s.add_mutually_exclusive_group()
    g.add_argument("--x", action="store_true")
    g.add_argument("-y", action="store_false")
    nested = sub.add_parser("tracks").add_subparsers(dest="track", required=False)
    nested.add_parser("gc").add_argument("--fasta", required=True)
    return p


def test_contract_of_toy_parser():
    c = build_contract(toy_parser())
    assert c["command_required"] is True and list(c["commands"]) == ["run", "tracks"]
    run = c["commands"]["run"]
    assert run["option_strings"] == ["--clock", "--count", "--help", "--level", "--mode", "--tag", "--x", "-h", "-n",
                                     "-y"]
    assert run["positionals"] == [{"dest": "infile", "flags": [], "required": True, "action": "store", "nargs": None,
                                   "type": None, "choices": None, "default": None}]
    o = run["options"]
    assert {k: v["action"] for k, v in o.items()} == {
        "--help": "help", "--count": "count", "--mode": "store_const", "--tag": "extend", "--level": "store",
        "--clock": "store", "--x": "store_true", "-y": "store_false"}
    assert o["--count"]["flags"] == ["--count", "-n"] and o["--count"]["default"] == 0
    assert o["--mode"]["const"] == "fast"
    assert o["--tag"]["nargs"] == "+" and o["--tag"]["default"] == ["a"]
    assert o["--level"]["type"] == "float" and o["--level"]["choices"] == [0.5, 1.0]
    assert "default" not in o["--clock"]
    assert run["mutually_exclusive"] == [{"required": False, "options": ["--x", "-y"]}]
    tracks = c["commands"]["tracks"]
    assert tracks["command_required"] is False and list(tracks["commands"]) == ["gc"]
    assert tracks["commands"]["gc"]["options"]["--fasta"]["required"] is True
    assert c["options"]["--version"]["version"] == "toy 1.2"  # no version given: nothing masked
    assert json.loads(render(c)) == c


def test_version_template_masks_the_version():
    """A changed --version format is a contract change (versions.yml parses it), a version bump is not."""
    v12 = build_contract(toy_parser(version="toy 1.2"), version="1.2")
    assert v12["options"]["--version"]["version"] == f"toy {VERSION_MASK}"
    assert contract_diff(v12, build_contract(toy_parser(version="toy 1.3"), version="1.3")) == []
    changed = build_contract(toy_parser(version="toy v1.2"), version="1.2")
    assert contract_diff(v12, changed) == [f'mycnv --version version: "toy {VERSION_MASK}" -> "toy v{VERSION_MASK}"']
    assert affected_modules(v12, changed) == [TOP_LEVEL_MODULES]


@pytest.mark.parametrize("new_default, shown", [(0.0, "0 -> 0.0"), (False, "0 -> false")])
def test_contract_diff_is_type_exact(new_default, shown):
    """0 == 0.0 == False in Python; in the snapshot (and for a Nextflow module) they are different defaults."""
    old, new = build_contract(toy_parser()), build_contract(toy_parser(count_default=new_default))
    assert old == new  # what a plain `==` comparison would have let through
    assert render(old) != render(new)
    assert contract_diff(old, new) == [f"mycnv run --count default: {shown}"]
    assert affected_modules(old, new) == ["`run`: no module calls it yet"]


def test_exact_flags_reject_a_renamed_flag():
    """A flag renamed to a longer one (--tag -> --tag-name) must fail a module call that still uses the old flag."""
    argv = ["run", "in.txt", "--tag", "b"]
    renamed = toy_parser(tag_flag="--tag-name")
    assert renamed.parse_args(argv).tag == ["a", "b"]  # argparse's prefix matching hides the rename
    with pytest.raises(SystemExit):
        exact_flags(toy_parser(tag_flag="--tag-name")).parse_args(argv)
    assert exact_flags(toy_parser()).parse_args(argv).tag == ["a", "b"]


def test_stale_install_is_named(tmp_path):
    """A test that imports another, different mycnv than src/ says so instead of 'rerun the script'.

    The workbench case: the env's editable install points at the kumodev checkout, the tests run in a worktree.
    """
    copy_ = tmp_path / "site-packages" / "mycnv"
    shutil.copytree(SRC_PKG, copy_, ignore=shutil.ignore_patterns("__pycache__"))
    installed = types.SimpleNamespace(__file__=str(copy_ / "__init__.py"), __version__="0.9")
    checkout = types.SimpleNamespace(__file__=str(SRC_PKG / "__init__.py"), __version__="1.0")
    assert not stale_install(checkout)
    assert not stale_install(installed)  # a non-editable install of this checkout (CI, Docker test stage)
    (copy_ / "cli.py").write_text((copy_ / "cli.py").read_text() + "\n# older\n")
    assert stale_install(installed)
    expected = build_contract()
    message = mismatch_message(expected, expected, module=installed)
    assert f"(mycnv 0.9 from {copy_.resolve()})" in message
    assert "stale install" in message and "PYTHONPATH=src" in message
    # PYTHONPATH first; a reinstall only in a private env: in the shared workbench env it repoints the editable install
    # of the kumodev checkout for every agent and worktree
    assert message.index("PYTHONPATH=src") < message.index("pip install -e .")
    assert "only in a private env, never in a shared one" in message
    assert "rerun `python tests/make_cli_contract.py`" not in message
    message = mismatch_message(expected, expected, module=checkout)
    assert f"(mycnv 1.0 from {SRC_PKG.resolve()})" in message
    assert "rerun `python tests/make_cli_contract.py`" in message and "stale" not in message


def test_mismatch_message_says_what_to_do():
    """The failure text: the changed entries, the rerun instruction and the modules to review (design section 3)."""
    old = build_contract()
    new = copy.deepcopy(old)
    new["commands"]["pon-coverage"]["options"]["--gc"]["required"] = False
    message = mismatch_message(old, new, module=types.SimpleNamespace(__file__=str(SRC_PKG / "__init__.py"),
                                                                      __version__=mycnv.__version__))
    assert "mycnv pon-coverage --gc required: true -> false" in message
    assert "rerun `python tests/make_cli_contract.py`" in message
    assert "review the kumo-variant-calling-nf modules that call the changed flags (design section 3" in message
    assert f"  {NF_MODULES['pon-coverage']}" in message


def test_contract_diff_names_the_change():
    old, new = build_contract(toy_parser()), build_contract(toy_parser(count_default=3))
    assert contract_diff(old, old) == []
    assert contract_diff(old, new) == ["mycnv run --count default: 0 -> 3"]
    del new["commands"]["run"]["options"]["--mode"]
    new["commands"]["run"]["option_strings"].remove("--mode")
    new["commands"]["run"]["options"]["--new"] = {"action": "store_true"}
    new["commands"]["run"]["option_strings"].append("--new")
    assert contract_diff(old, new) == [
        "mycnv run option_strings: added --new; removed --mode",
        "mycnv run --count default: 0 -> 3",
        "mycnv run --mode: removed",
        'mycnv run --new: added {"action": "store_true"}',
    ]
    assert affected_modules(old, new) == ["`run`: no module calls it yet"]
