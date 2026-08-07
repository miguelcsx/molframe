# Space-group catalogue

`space-groups.bin` is an offline catalogue of the 530 Hall settings covering
all 230 crystallographic space-group types. It stores exact integer rotations,
translations as twelfths, International Tables numbers, Hermann–Mauguin names,
Hall symbols and setting choices.

The data is generated from spglib 2.7.0 with:

```bash
python scripts/generate-space-groups.py \
  crates/pdbiox-xtal/data/space-groups.bin
```

The binary format and decoder are owned by pdbiox. The source records come from
spglib's database and retain its BSD-3-Clause notice in `LICENSE.spglib`.
