# Space-group catalogue

`space-groups.bin` is an offline catalogue of the 530 Hall settings covering
all 230 crystallographic space-group types. It stores exact integer rotations,
translations as twelfths, International Tables numbers, Hermann–Mauguin names,
Hall symbols and setting choices.

The records were derived from spglib 2.7.0 and are committed as a payload — the
file is data, not a build step, and there is no script to regenerate it.

The binary format and decoder are owned by molframe: an 8-byte `PDBXSG01` magic,
then 530 length-prefixed setting records. The source records come from spglib's
database and retain its BSD-3-Clause notice in `LICENSE.spglib`.
