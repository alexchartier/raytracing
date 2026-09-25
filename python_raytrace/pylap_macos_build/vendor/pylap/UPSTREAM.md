# Vendored PyLap raytrace source

Source: HamSCI/PyLap, commit `4b32d0a5257179699ecf3ea72fb4d49c8b9aefde`
(https://github.com/HamSCI/PyLap).
License: MIT; see [LICENSE](LICENSE).

This directory contains `modules/source/raytrace_3d.c`, its six
`modules/source/common/*.c` helpers, and `modules/include/pharlap.h`.
Generated binaries, examples, and unrelated PyLap modules are excluded.

Local changes to `raytrace_3d.c`:

- Parse the 10-argument cached-grid form with the ray state vector in argument
  10, matching the MATLAB API.
- Return an error when grid validation or loading fails, and correct the two
  helper return types for current macOS Clang.

The PHaRLAP libraries and headers are separate dependencies supplied by
`PHARLAP_HOME`.
