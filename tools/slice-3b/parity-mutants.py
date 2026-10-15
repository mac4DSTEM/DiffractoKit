#!/usr/bin/env python3
"""Sensitivity of the bit-identity gate (registration section 5, gate 3): each EllipseCalibration mutant of mutants.py is compiled
(-O, whole module, the same flags as run.sh's kit-check) with the apps-parity main.swift and its dump must DIFFER from the apps' -O
dump (.build/apps-parity-3b/dump-O.json). A mutant whose dump is byte-identical is not seen by the 35 rings.

  python3 tools/slice-3b/parity-mutants.py <scratch-dir>
"""
import importlib.util
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("mutants", ROOT / "tools/slice-3b/mutants.py")
mutants = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mutants)

scratch = pathlib.Path(sys.argv[1])
scratch.mkdir(parents=True, exist_ok=True)
work = ROOT / ".build/apps-parity-3b"
reference = (work / "dump-O.json" if (work / "dump-O.json").exists() else work / "dump.json").read_bytes()
geo = ROOT / "Sources/DiffractoGeometry"
procs = []
for ident, rel, old, new, _flt, what in mutants.MUTANTS:
    if not rel.endswith("EllipseCalibration.swift"):
        continue
    d = scratch / ident
    d.mkdir(exist_ok=True)
    text = (ROOT / rel).read_text()
    assert text.count(old) == 1, ident
    (d / "EllipseCalibration.swift").write_text(text.replace(old, new))
    files = [str(d / "EllipseCalibration.swift")] + [str(geo / f) for f in ("DiffractionPattern.swift", "RingCalibration.swift", "ElectronWavelength.swift")]
    files += [str(ROOT / "tools/slice-3b/apps-parity/kit-shim.swift"), str(ROOT / "tools/slice-3b/apps-parity/main.swift")]
    subprocess.run(["swiftc", "-O", "-wmo", "-swift-version", "6", "-package-name", "x", "-o", str(d / "tool")] + files, check=True)
    procs.append((ident, what, d, subprocess.Popen([str(d / "tool"), str(d / "dump.json")], stderr=subprocess.DEVNULL)))
seen = 0
for ident, what, d, p in procs:
    p.wait()
    differs = (d / "dump.json").read_bytes() != reference
    seen += differs
    print("%s [%s]: %s" % (ident, what, "dump DIFFERS (seen by the parity gate)" if differs else "dump identical (NOT seen)"), flush=True)
print("parity gate: %d of %d EllipseCalibration mutants change the dump" % (seen, len(procs)))
