#!/usr/bin/env python3
"""Build the bright-star catalog files from the HYG database.

    python3 scripts/build_catalog.py --source source/hyg_v44.csv.gz --out data/

Writes data/stars-full.bin, data/stars-balanced.bin and data/manifest.json. Each .bin is a flat array of
little-endian float32, 5 values per star (x, y, z, mag, bv), 20 bytes per record, no header; see README.md.

How a cut is made: sort every star except the Sun by magnitude (then by HYG id), take the brightest N (the target
count), then keep every further star exactly as bright as the Nth. So each file is "all stars at or brighter than
its magnitude limit", and the balanced file is a prefix of the full one.

Standard library only. The same source file gives byte-identical outputs: there are no timestamps, hostnames or
paths in the manifest.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import struct
import sys
from pathlib import Path

# The pinned source. A new HYG release is a reviewed change to these values, not an automatic upgrade.
SOURCE = {
    'catalog': 'HYG',
    'version': '4.4',
    'repository': 'https://codeberg.org/astronexus/hyg',
    'commit': '53e3df311869e813ace5f1ad2ec4ce909f13256c',
    'path': 'data/hyg/CURRENT/hyg_v44.csv.gz',
    'url': 'https://codeberg.org/astronexus/hyg/media/commit/53e3df311869e813ace5f1ad2ec4ce909f13256c/data/hyg/CURRENT/hyg_v44.csv.gz',
    'sha256': '00b349893b9a53106dd488d8371e8d2fa586043e500bb3cdb8bff3931682197d',
    'bytes': 13636362,
    'license': 'CC-BY-SA-4.0',
}
TARGETS = {'full': 6000, 'balanced': 1000}
FIELDS = ('x', 'y', 'z', 'mag', 'bv')
RECORD = struct.Struct('<5f')
STRIDE = RECORD.size  # 20
MAX_SOURCE_BYTES = 100_000_000
NUMBER = re.compile(r'^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$')


class BuildError(Exception):
    pass


def _number(text: str, field: str, row: int) -> float:
    t = text.strip()
    if not NUMBER.match(t):
        raise BuildError(f'row {row}: {field} is not a number: {t!r}')
    v = float(t)
    if not math.isfinite(v) or not math.isfinite(struct.unpack('<f', struct.pack('<f', v))[0]):
        raise BuildError(f'row {row}: {field} is out of float32 range')
    return v


def read_source(path: Path) -> bytes:
    data = path.read_bytes()
    if len(data) > MAX_SOURCE_BYTES:
        raise BuildError(f'source is larger than {MAX_SOURCE_BYTES} bytes')
    return data


def parse_stars(raw: bytes) -> list[tuple[float, int, float, float, float]]:
    """(mag, id, ra_hours, dec_degrees, bv) for every star but the Sun, brightest first, then by HYG id."""
    text = (gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw).decode('utf-8-sig')
    reader = csv.reader(io.StringIO(text, newline=''))
    header = next(reader, None)
    if not header:
        raise BuildError('source has no header')
    col = {name: i for i, name in enumerate(header)}
    for name in ('id', 'proper', 'ra', 'dec', 'mag', 'ci'):
        if name not in col:
            raise BuildError(f'source has no {name} column')
    stars, ids, suns = [], set(), 0
    for n, r in enumerate(reader, start=2):
        if not r:
            continue
        if len(r) != len(header):
            raise BuildError(f'row {n}: {len(r)} fields, header has {len(header)}')
        sid = _number(r[col['id']], 'id', n)
        if sid != int(sid) or sid < 0 or sid in ids:
            raise BuildError(f'row {n}: bad or duplicate id')
        sid = int(sid)
        ids.add(sid)
        proper = r[col['proper']].strip()
        if sid == 0:
            if proper != 'Sol':
                raise BuildError('id 0 is not the Sun')
            suns += 1
            continue
        if proper.lower() == 'sol':
            raise BuildError(f'row {n}: a second Sun')
        ra, dec, mag = (_number(r[col[f]], f, n) for f in ('ra', 'dec', 'mag'))
        if not (0 <= ra < 24 and -90 <= dec <= 90):
            raise BuildError(f'row {n}: position out of range')
        ci = r[col['ci']]
        bv = math.nan if ci.strip() == '' else _number(ci, 'ci', n)
        stars.append((mag, sid, ra, dec, bv))
    if suns != 1:
        raise BuildError('source has no Sun record (id 0)')
    stars.sort(key=lambda s: (s[0], s[1]))
    return stars


def cut(stars, target: int):
    """The brightest `target` stars plus any further stars as bright as the last one; and that magnitude."""
    if not 0 < target <= len(stars):
        raise BuildError(f'target {target} needs at least that many stars (have {len(stars)})')
    limit = stars[target - 1][0]
    end = target
    while end < len(stars) and stars[end][0] == limit:
        end += 1
    return limit, stars[:end]


def pack(stars) -> bytes:
    out = bytearray()
    for mag, _sid, ra, dec, bv in stars:
        a, d = ra * math.pi / 12, dec * math.pi / 180
        c = math.cos(d)
        out += RECORD.pack(c * math.cos(a), c * math.sin(a), math.sin(d), mag, bv)
    return bytes(out)


def build(source: Path, out: Path, expected_sha256: str | None = SOURCE['sha256'], targets: dict[str, int] = TARGETS,
          source_meta: dict | None = None) -> dict:
    raw = read_source(source)
    got = hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and got != expected_sha256:
        raise BuildError(f'source sha256 {got} is not the pinned {expected_sha256}')
    stars = parse_stars(raw)
    if targets['full'] <= targets['balanced']:
        raise BuildError('the full target must be larger than the balanced one')
    files, blobs = {}, {}
    for name in ('full', 'balanced'):
        limit, subset = cut(stars, targets[name])
        blob = pack(subset)
        blobs[name] = blob
        files[name] = {
            'file': f'stars-{name}.bin',
            'targetCount': targets[name],
            'count': len(subset),
            'magnitudeLimit': limit,
            'brightestMagnitude': subset[0][0],
            'missingBv': sum(1 for s in subset if math.isnan(s[4])),
            'bytes': len(blob),
            'sha256': hashlib.sha256(blob).hexdigest(),
        }
    meta = dict(source_meta if source_meta is not None else SOURCE)
    meta['sha256'] = got
    meta['bytes'] = len(raw)
    manifest = {
        'manifestVersion': 1,
        'source': meta,
        'license': {'data': 'CC-BY-SA-4.0', 'code': 'MIT'},
        'epoch': 'J2000',
        'frame': 'equatorial: +x toward RA 0h Dec 0, +y toward RA 6h Dec 0, +z toward the north celestial pole',
        'record': {'fields': list(FIELDS), 'type': 'float32', 'byteOrder': 'little-endian', 'stride': STRIDE,
                   'header': None, 'missingBv': 'NaN'},
        'order': 'magnitude ascending, then HYG id ascending',
        'cut': 'the brightest targetCount stars, then every star as bright as the last: magnitude <= magnitudeLimit',
        'excluded': ['the Sun (HYG id 0)'],
        'eligibleStars': len(stars),
        'files': files,
        'recipe': 'python3 scripts/build_catalog.py --source <hyg csv.gz> --out data/',
    }
    out.mkdir(parents=True, exist_ok=True)
    for name, blob in blobs.items():
        (out / files[name]['file']).write_bytes(blob)
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    return manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description='Build the bright-star catalog files from the HYG database.')
    ap.add_argument('--source', type=Path, required=True, help='the HYG csv or csv.gz (pinned in SOURCE)')
    ap.add_argument('--out', type=Path, default=Path('data'))
    a = ap.parse_args(argv)
    try:
        m = build(a.source, a.out)
    except (BuildError, OSError, UnicodeDecodeError, csv.Error) as e:
        print(f'build failed: {e}', file=sys.stderr)
        return 1
    for name, f in m['files'].items():
        print(f'{name}: {f["count"]} stars, magnitude <= {f["magnitudeLimit"]}, {f["bytes"]} bytes, sha256 {f["sha256"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
