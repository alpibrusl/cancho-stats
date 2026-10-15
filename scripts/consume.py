#!/usr/bin/env python3
"""G9 (docs/numerics.md section 7.4): the consumption proof. A scratch package OUTSIDE this repository with every `[dependencies.*]` block of cancho.toml, `cancho install`, and a program that
imports every pinned module and calls one function of each family, comparing the values with numbers written here (from the definitions, not from the code under test).

    python3 scripts/consume.py [--cancho PATH] [--mutants]

Exit 0: the install and the build succeed and every value is the expected one. Exit 1: an install, build or value failure. Exit 2: the harness could not run.
`--mutants`: a block removed, a wrong expected value, and one byte of a store changed in a scratch clone (cancho install must refuse it): each must make the proof fail.
"""
import hashlib
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))

MAIN = r'''edition 5;

import std.io;
import std.math;
import ml.rng;
import ml.exact;
import ml.fbits;
import ml.fsq;
import ml.limits;
import ml.refuse;
import ml.matrix;
import linalg.order;
import linalg.tensor;
import linalg.dot;
import linalg.axpy;
import linalg.gemv;
import linalg.syrk;
import linalg.gemm;
import linalg.chol;
import linalg.par;
import metrics.rules;
import metrics.ratio;
import metrics.classify;
import metrics.sums;
import metrics.regress;
import metrics.rank;
import metrics.splits;
import facc;
import flt;
import dec;
import pow5;
import reader;
import toolbox.cli;
import toolbox.describe;
import toolbox.fail;
import toolbox.limit;
import toolbox.lines;
import toolbox.out;
import toolbox.path;
import toolbox.place;
import toolbox.text;

fn hex16[&i](io: &!i Io, bits: int) -> [io_write] int {
    var k = 0;
    while k < 16 {
        let nib = bits >> 60 - 4 * k & 15;
        if nib < 10 {
            putchar(io, '0' + nib);
        } else {
            putchar(io, 'a' + nib - 10);
        }
        k = k + 1;
    }
    return 0;
}

fn line[&i](io: &!i Io, name: &static [byte], bits: int) -> [io_write] int {
    io.write_all(io, name);
    io.space(io);
    hex16(io, bits);
    io.newline(io);
    return 0;
}

fn main(world: World) -> [] int {
    let Split { io, ffi, fs, heap, args, net, clock } = split(world);
    release(ffi);
    release(net);
    release(fs);
    release(clock);
    release(args);
    var nine = 0.0;
    borrow mut heap as &!h in {
        let wide = box_slice(h, fsq.size() + 1, 0);
        borrow mut wide as &!w in {
            fsq.add_square(contents(w), bits_of(3.0), 0);
            let (v, st) = fsq.finish(contents(w));
            nine = v;
        }
        unbox_slice(h, wide);
    }
    release(heap);
    borrow mut io as &!w in {
        line(w, "rng.first", rng.output(rng.advance(0)));
        line(w, "rng.derive", rng.derive(rng.derive(42, 1), 0));
        line(w, "exact.pow10", exact.pow10(5));
        line(w, "fbits.roundtrip", bits_of(fbits.decode(bits_of(1.5))));
        line(w, "fsq.square", bits_of(nine));
        line(w, "refuse.count", refuse.count());
        line(w, "ratio", bits_of(ratio.ratio(1, 4)));
        line(w, "facc.width", facc.width());
        line(w, "flt.roundtrip", bits_of(flt.float_of_key(flt.key_of(1.5))));
        line(w, "order.lanes", order.lanes());
        line(w, "tensor.rmax", tensor.rmax());
        line(w, "chol.pivot_tol", bits_of(chol.pivot_tol()));
        line(w, "rules.count", rules.count());
    }
    release(io);
    return 0;
}
'''

EXPECT = {
    "rng.first": 0xE220A8397B1DCDAF,  # SplitMix64 from seed 0: the first output (the C reference, cancho-ml docs/rng.md)
    "rng.derive": 0xD8005CDD08A0D146,  # derive(derive(42, 1), 0), same document
    "exact.pow10": 100000,
    "fbits.roundtrip": 0x3FF8000000000000,  # 1.5
    "fsq.square": 0x4022000000000000,  # 9.0
    "ratio": 0x3FD0000000000000,  # 0.25
    "facc.width": 73,  # 72 limbs and a counter (cancho-table docs/numbers.md)
    "flt.roundtrip": 0x3FF8000000000000,
}
# only run, not compared with a number written here: the value is the project's own constant (listed so the proof reaches the module)
REACHED = ["refuse.count", "order.lanes", "tensor.rmax", "chol.pivot_tol", "rules.count"]


def run(cmd, cwd, **kw):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, **kw)


def proof(toml_text, expect, label):
    """-> (ok, message, sha256 of the output)"""
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        text = re.sub(r"\n\[\[bin\]\].*", "", toml_text, flags=re.S)
        text += '\n[[bin]]\nname = "consume"\nsources = ["src"]\nstd = true\nout = "build/consume"\n'
        (tmp / "cancho.toml").write_text(text)
        (tmp / "src").mkdir()
        (tmp / "src" / "main.cho").write_text(MAIN)
        r = run([CANCHO, "install"], tmp)
        if r.returncode != 0:
            return False, "%s: cancho install failed: %s" % (label, (r.stderr or r.stdout)[-300:]), ""
        r = run([CANCHO, "build", "--bin", "consume"], tmp)
        if r.returncode != 0:
            return False, "%s: build failed: %s" % (label, (r.stderr or r.stdout)[-400:]), ""
        r = run([str(tmp / "build" / "consume")], tmp)
        if r.returncode != 0:
            return False, "%s: the program exited %d" % (label, r.returncode), ""
        got = {}
        for ln in r.stdout.strip().split("\n"):
            k, v = ln.split()
            got[k] = int(v, 16)
        for k, want in expect.items():
            if got.get(k) != want:
                return False, "%s: %s is %s, expected %s" % (label, k, hex(got.get(k, -1)), hex(want)), ""
        for k in REACHED:
            if k not in got:
                return False, "%s: %s was not reached" % (label, k), ""
        return True, "%s: %d values equal the expected, %d modules reached" % (label, len(expect), len(REACHED)), hashlib.sha256(r.stdout.encode()).hexdigest()


def main():
    global CANCHO
    if "--cancho" in sys.argv:
        CANCHO = sys.argv[sys.argv.index("--cancho") + 1]
    if not pathlib.Path(CANCHO).exists():
        print("consume: no compiler at %s" % CANCHO)
        return 2
    toml = (ROOT / "cancho.toml").read_text()
    ok, msg, digest = proof(toml, EXPECT, "proof")
    print(msg + ("  sha256 " + digest if digest else ""))
    if not ok:
        return 1
    if "--mutants" not in sys.argv:
        return 0
    survived = []
    # 1. a block removed (facc: no other store requires it, whereas the `requires` of a numeric store of cancho-ml would bring ml.rng in anyway)
    m1 = re.sub(r"\[dependencies\.facc\]\n(?:[^\n\[]*\n)*", "", toml, count=1)
    ok1, msg1, _ = proof(m1, EXPECT, "mutant block-removed")
    print(msg1 if not ok1 else "mutant block-removed SURVIVED")
    if ok1:
        survived.append("block-removed")
    # 2. a wrong expected value
    bad = dict(EXPECT)
    bad["fsq.square"] = 0x4022000000000001
    ok2, msg2, _ = proof(toml, bad, "mutant wrong-value")
    print(msg2 if not ok2 else "mutant wrong-value SURVIVED")
    if ok2:
        survived.append("wrong-value")
    # 3. one byte of a store changed in a scratch clone of the numeric stores, pinned at the new commit
    m = re.search(r'git = "(https://github.com/alpibrusl/cancho-ml)"', toml)
    if m:
        src = m.group(1)
        with tempfile.TemporaryDirectory() as t:
            clone = pathlib.Path(t) / "stores"
            run(["git", "clone", "-q", src, str(clone)], t)
            target = next((clone / ".cancho-vcs" / "ml.rng").rglob("*"), None)
            files = [p for p in (clone / ".cancho-vcs" / "ml.rng").rglob("*") if p.is_file() and p.stat().st_size > 20]
            p = files[0]
            data = bytearray(p.read_bytes())
            data[len(data) // 2] ^= 1
            p.write_bytes(bytes(data))
            run(["git", "-c", "user.name=m", "-c", "user.email=m@m", "commit", "-q", "-a", "-m", "tamper"], clone)
            rev = run(["git", "rev-parse", "HEAD"], clone).stdout.strip()
            old = re.search(r'rev = "([0-9a-f]{40})"\npath = "\.cancho-vcs/ml\.rng"', toml).group(1)
            m3 = toml.replace(src, "file://" + str(clone)).replace(old, rev)
            ok3, msg3, _ = proof(m3, EXPECT, "mutant store-byte-changed")
            print(msg3 if not ok3 else "mutant store-byte-changed SURVIVED")
            if ok3:
                survived.append("store-byte-changed")
    print("consume --mutants: %d surviving" % len(survived))
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())
