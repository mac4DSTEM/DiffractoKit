# Security policy

## Supported versions

DiffractoKit is at version 0.x. Security fixes go into the latest release only.

## Reporting a vulnerability

Please do not open a public issue for a security problem. Use GitHub's private vulnerability reporting for this repository
(the "Report a vulnerability" button under the Security tab). Include the version, a description, and if possible a small
input (for example a CIF file) that reproduces it. You can expect an acknowledgement within a few days.

## Scope

The library parses CIF text and reads numeric image data supplied by the caller. Reports about crashes, traps, unbounded
memory or time on malformed input are in scope.
