# readme-consumer

A separate SwiftPM package that depends on DiffractoKit **by path** (`.package(path: "../..")`) and compiles, from outside
the package, exactly the program in the main README's "First ten minutes". Its target, `Sources/readme-consumer`, is a symbolic
link to the kit's `Examples/diffracto-example`, so the program exists once: `Examples/diffracto-example/main.swift` is that code
and nothing else. `run.sh` cuts the program and its stated output out of `README.md` and fails unless

1. the example's `main.swift` equals the README's program, byte for byte;
2. the package builds in Debug with no warning (and no test-only target is built);
3. the program's output equals the README's output block, line for line.

When the README's code changes, copy it to `Examples/diffracto-example/main.swift` and the new output
into the README's output block. [`tools/consumer-check`](../consumer-check/README.md) is the wider guard of the public API.

```
tools/readme-consumer/run.sh
READMEC_JOBS=2 tools/readme-consumer/run.sh      # swift build -j 2
```
