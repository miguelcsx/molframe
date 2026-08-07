#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Miguel Cárdenas
# SPDX-License-Identifier: MIT OR Apache-2.0
#
# Runs `reuse lint` (https://reuse.software) against the working tree. Exits
# zero on compliance, non-zero on any REUSE finding. verify.sh guards the
# presence of the tool before calling this script.
set -uo pipefail

cd "$(dirname "$0")/.."

reuse lint
