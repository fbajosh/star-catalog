"""Tests for build_catalog.py. Standard library only; no download needed.

    python3 -m unittest discover scripts

Synthetic HYG-format rows test the build. If data/ holds a built catalog, it is also checked against its manifest.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import struct
import tempfile
import unittest
from pathlib import Path

import build_catalog as bc

HEADER = ('id,hip,hd,hr,gl,bf,proper,ra,dec,dist,pmra,pmdec,rv,mag,absmag,spect,ci,x,y,z,vx,vy,vz,rarad,decrad,'
          'pmrarad,pmdecrad,bayer,flam,con,comp,comp_primary,base,lum,var,var_min,var_max')
NCOLS = len(HEADER.split(','))
DATA = Path(__file__).resolve().parent.parent / 'data'


def row(sid, ra, dec, mag, ci='0.5', proper=''):
    cells = [''] * NCOLS
    cells[0], cells[6], cells[7], cells[8], cells[13], cells[16] = str(sid), proper, str(ra), str(dec), str(mag), ci
    return ','.join(cells)


def fixture(rows, gz=True):
    text = '\n'.join([HEADER, row(0, 0, 0, -26.7, '0.656', 'Sol'), *rows]) + '\n'
    raw = text.encode()
    return gzip.compress(raw, mtime=0) if gz else raw


# Axis stars, ties at magnitude 3.0, and fillers. Targets in the tests are small (full 6, balanced 2).
BASE = [
    row(1, 0, 0, 1.0, '-0.2'),        # +x
    row(2, 6, 0, 2.0, '0.1'),         # +y
    row(3, 12, 89.999999, 2.5, ''),   # toward +z; blank B-V
    row(4, 18, -45, 3.0, '1.0'),
    row(5, 3, 30, 3.0, '0.4'),        # tie with 4
    row(6, 9, -60, 3.0, '0.2'),       # tie with 4 and 5
    row(7, 21, 10, 4.0, '0.3'),
    row(8, 1, 1, 5.0, '0.9'),
    row(9, 2, 2, 6.0, '1.2'),
]


class Build(unittest.TestCase):
    def build(self, rows=BASE, targets=None, gz=True):
        tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        src = tmp / ('in.csv.gz' if gz else 'in.csv')
        src.write_bytes(fixture(rows, gz))
        out = tmp / 'out'
        m = bc.build(src, out, expected_sha256=None, targets=targets or {'full': 5, 'balanced': 2},
                     source_meta={'catalog': 'fixture', 'url': 'https://example.invalid/fixture'})
        return m, out

    def records(self, path):
        b = path.read_bytes()
        self.assertEqual(len(b) % 20, 0)
        return [struct.unpack_from('<5f', b, i) for i in range(0, len(b), 20)]

    def test_layout_axes_and_unit_vectors(self):
        m, out = self.build()
        recs = self.records(out / 'stars-full.bin')
        self.assertEqual(m['record']['stride'], 20)
        self.assertEqual(m['record']['fields'], ['x', 'y', 'z', 'mag', 'bv'])
        x, y, z, mag, bv = recs[0]
        self.assertAlmostEqual(x, 1, places=6); self.assertAlmostEqual(y, 0, places=6); self.assertAlmostEqual(z, 0, places=6)
        self.assertEqual((mag, bv), (1.0, struct.unpack('<f', struct.pack('<f', -0.2))[0]))
        self.assertAlmostEqual(recs[1][1], 1, places=6)   # RA 6h is +y
        self.assertAlmostEqual(recs[2][2], 1, places=6)   # Dec 90 is +z
        self.assertTrue(math.isnan(recs[2][4]))           # blank B-V
        for r in recs:
            self.assertAlmostEqual(math.hypot(r[0], r[1], r[2]), 1, delta=1e-6)

    def test_cut_keeps_ties_and_reports_the_limit(self):
        m, out = self.build(targets={'full': 4, 'balanced': 2})
        f = m['files']['full']
        # The 4th brightest is magnitude 3.0, shared by ids 4, 5 and 6: all three are kept.
        self.assertEqual((f['targetCount'], f['count'], f['magnitudeLimit']), (4, 6, 3.0))
        mags = [r[3] for r in self.records(out / 'stars-full.bin')]
        self.assertEqual(mags, sorted(mags))
        self.assertEqual(m['files']['balanced']['magnitudeLimit'], 2.0)
        self.assertEqual(m['files']['full']['missingBv'], 1)

    def test_ties_break_by_id_and_balanced_is_a_prefix(self):
        _, out = self.build(targets={'full': 6, 'balanced': 2})
        full, bal = (out / 'stars-full.bin').read_bytes(), (out / 'stars-balanced.bin').read_bytes()
        self.assertEqual(full[:len(bal)], bal)
        # ids 4, 5, 6 (RA 18h, 3h, 9h) in id order.
        ras = [math.degrees(math.atan2(r[1], r[0])) % 360 for r in self.records(out / 'stars-full.bin')[3:6]]
        self.assertEqual([round(a) for a in ras], [270, 45, 135])

    def test_the_sun_is_left_out(self):
        _, out = self.build(targets={'full': 9, 'balanced': 1})
        mags = [r[3] for r in self.records(out / 'stars-full.bin')]
        self.assertEqual(len(mags), 9)
        self.assertGreater(min(mags), -2)

    def test_same_source_same_bytes(self):
        _, a = self.build()
        _, b = self.build()
        for name in ('stars-full.bin', 'stars-balanced.bin', 'manifest.json'):
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)

    def test_plain_csv_matches_gzip(self):
        _, a = self.build(gz=True)
        _, b = self.build(gz=False)
        self.assertEqual((a / 'stars-full.bin').read_bytes(), (b / 'stars-full.bin').read_bytes())

    def test_manifest_hashes_and_no_local_paths(self):
        m, out = self.build()
        for f in m['files'].values():
            self.assertEqual(hashlib.sha256((out / f['file']).read_bytes()).hexdigest(), f['sha256'])
        text = (out / 'manifest.json').read_text()
        self.assertNotIn(str(out.parent), text)
        self.assertNotIn(str(Path.home()), text)

    def test_bad_input_fails(self):
        cases = {
            'non-numeric B-V': BASE + [row(10, 1, 1, 7.0, 'x')],
            'duplicate id': BASE + [row(1, 1, 1, 7.0)],
            'RA out of range': BASE + [row(10, 24, 0, 7.0)],
            'second Sun': BASE + [row(10, 1, 1, 7.0, proper='Sol')],
            'short row': BASE + ['10,,,'],
        }
        for name, rows in cases.items():
            with self.subTest(name), self.assertRaises(bc.BuildError):
                self.build(rows)
        with self.assertRaises(bc.BuildError):
            self.build(targets={'full': 20, 'balanced': 2})   # more than there are

    def test_pinned_hash_is_enforced(self):
        tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (tmp / 'in.csv.gz').write_bytes(fixture(BASE))
        with self.assertRaises(bc.BuildError):
            bc.build(tmp / 'in.csv.gz', tmp / 'out')
        self.assertFalse((tmp / 'out').exists())


@unittest.skipUnless((DATA / 'manifest.json').exists(), 'no built catalog in data/')
class CommittedData(unittest.TestCase):
    def test_files_match_the_manifest(self):
        m = json.loads((DATA / 'manifest.json').read_text())
        self.assertEqual(m['source']['sha256'], bc.SOURCE['sha256'])
        prev = None
        for name in ('full', 'balanced'):
            f = m['files'][name]
            b = (DATA / f['file']).read_bytes()
            self.assertEqual((len(b), len(b) // 20, hashlib.sha256(b).hexdigest()), (f['bytes'], f['count'], f['sha256']))
            mags = [struct.unpack_from('<f', b, i + 12)[0] for i in range(0, len(b), 20)]
            self.assertEqual(mags, sorted(mags))
            self.assertEqual(mags[-1], struct.unpack('<f', struct.pack('<f', f['magnitudeLimit']))[0])
            if prev is not None:
                self.assertEqual(prev[:len(b)], b)   # balanced is a prefix of full
            prev = b

    def test_sirius_is_first(self):
        # Sirius: J2000 RA 06h45m08.9s, Dec -16°42'58" (101.2872°, -16.7161°); HYG gives V -1.44, B-V 0.009.
        x, y, z, mag, bv = struct.unpack_from('<5f', (DATA / 'stars-balanced.bin').read_bytes(), 0)
        self.assertAlmostEqual(math.degrees(math.atan2(y, x)) % 360, 101.2872, delta=0.001)
        self.assertAlmostEqual(math.degrees(math.asin(z)), -16.7161, delta=0.001)
        self.assertAlmostEqual(mag, -1.44, places=5)
        self.assertAlmostEqual(bv, 0.009, places=5)


if __name__ == '__main__':
    unittest.main()
