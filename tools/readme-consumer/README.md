# readme-consumer

A separate SwiftPM package that depends on DiffractoKit **by path** (`.package(path: "../..")`) and compiles, from outside
the package, exactly the program in the main README's "First ten minutes". `Sources/readme-consumer/main.swift` is that code
and nothing else; so is `Examples/diffracto-example/main.swift`. `run.sh` cuts the program and its stated output out of
`README.md` and fails unless

1. both `main.swift` files equal the README's program, byte for byte;
2. the package builds in Debug with no warning (and no test-only target is built);
3. the program's output equals the README's output block, line for line.

When the README's code changes, copy it to both `main.swift` files (the example's and this package's) and the new output
into the README's output block. [`tools/consumer-check`](../consumer-check/README.md) is the wider guard of the public API.

```
tools/readme-consumer/run.sh
READMEC_JOBS=2 tools/readme-consumer/run.sh      # swift build -j 2
```
