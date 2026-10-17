# star-catalog

A compact bright-star catalog for drawing a realistic night sky in a WebGL globe. Each star is a fixed-size binary record (a direction, its brightness and its color). A renderer can load the file as one typed array without parsing.

The data is derived from the [HYG database](https://codeberg.org/astronexus/hyg) (CC BY-SA 4.0). This repository publishes the derived subset and the script that builds it, so anyone can rebuild or check it.

## Files

| Path | What it is |
|---|---|
| `data/stars-full.bin` | About 6,000 brightest stars (the magnitude limit is in `data/manifest.json`) |
| `data/stars-balanced.bin` | About 1,000 brightest stars, for lighter rendering |
| `data/manifest.json` | Source version and checksum, the cut for each file, star counts, sha256 of each output |
| `scripts/` | The conversion script and its test |

## Record format

Each file is a flat array of little-endian `float32`, **5 values per star, 20 bytes per record**, with no header:

| Offset | Field | Meaning |
|---|---|---|
| 0, 4, 8 | `x, y, z` | Unit vector toward the star, equatorial J2000 (x toward RA 0h Dec 0°, z toward the north celestial pole) |
| 12 | `mag` | Apparent visual magnitude (smaller is brighter) |
| 16 | `bv` | B−V color index (blue is negative, red is about 2); `NaN` where the source has none |

Stars are sorted brightest first, so a prefix of any file is also a valid, brighter cut.

The Sun is excluded. Proper motion is not applied: positions are at the source's epoch (J2000), which is accurate well below a pixel for this use.

## Rebuild

```sh
python3 scripts/build_catalog.py --source <hyg csv.gz> --out data/
python3 -m unittest discover scripts
```

The build is deterministic: the same source file gives byte-identical outputs, and `data/manifest.json` records the source's sha256.

## License

- **Data** (`data/`): [CC BY-SA 4.0](LICENSE-DATA), as a derivative of the HYG database. See [NOTICE](NOTICE) for attribution and the changes made.
- **Code** (`scripts/`): [MIT](LICENSE).
