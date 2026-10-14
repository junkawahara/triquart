#!/usr/bin/env python3
"""Compare fable5_1 (soda13 CLI) and astra6 (astra_cli) on many inputs.

Every case writes two Newick files and checks that both programs print the
same A..E counts for the triplet and the quartet distance.  Run through
run.sh, which builds both programs first.

  crosscheck.py FABLE_BIN ASTRA_BIN [--quick]
"""
import os, random, re, subprocess, sys, tempfile, threading

sys.setrecursionlimit(1_000_000)
threading.stack_size(512 * 1024 * 1024)

COUNTS = re.compile(r'A=(\d+) B=(\d+) C=(\d+) D=(\d+) E=(\d+)')


def ser(t, unary):
    """Serialize a nested-list tree; a 1-tuple is a unary node."""
    if isinstance(t, str):
        return t
    if isinstance(t, tuple):
        inner = ser(t[0], unary)
        return '(' + inner + ')' if unary else inner
    return '(' + ','.join(ser(c, unary) for c in t) + ')'


def decorate(s, rng):
    """Add random branch lengths."""
    return re.sub(r'([\w\)])(?=[,\)])',
                  lambda m: m.group(1) + (':%.2f' % rng.random() if rng.random() < .2 else ''), s)


def rtree(labels, shape, rng):
    nodes = labels[:]
    rng.shuffle(nodes)
    if shape == 'cat':
        cur = nodes[0]
        for x in nodes[1:]:
            cur = [cur, x]
        return cur
    if shape == 'star':
        return nodes
    maxd = {'bin': 2, 'mix': 5, 'high': 12}[shape]
    while len(nodes) > 1:
        k = min(len(nodes), rng.randint(2, maxd))
        rng.shuffle(nodes)
        grp, nodes = nodes[:k], nodes[k:]
        nodes.append((grp,) if rng.random() < 0.15 else grp)
    return nodes[0]


class Checker:
    def __init__(self, fable, astra, tmp):
        self.fable, self.astra, self.tmp = fable, astra, tmp
        self.total = self.bad = 0

    def path(self, name, text):
        p = os.path.join(self.tmp, name)
        with open(p, 'w') as f:
            f.write(text)
        return p

    def run(self, exe, *args):
        r = subprocess.run([exe, *args], capture_output=True, text=True)
        return r.returncode, COUNTS.findall(r.stdout)

    def compare(self, desc, fable_trees, astra_trees):
        """fable_trees may contain unary nodes; astra_trees must not."""
        self.total += 1
        f1, f2 = (self.path('f1.nwk', fable_trees[0]), self.path('f2.nwk', fable_trees[1]))
        a1, a2 = (self.path('a1.nwk', astra_trees[0]), self.path('a2.nwk', astra_trees[1]))
        frc, fc = self.run(self.fable, f1, f2)
        arc, ac = self.run(self.astra, a1, a2)
        if frc != 0 or arc != 0 or len(fc) != 2 or fc != ac:
            self.bad += 1
            print('MISMATCH', desc, 'fable:', frc, fc, 'astra:', arc, ac)
            print('  T1:', fable_trees[0][:300])
            print('  T2:', fable_trees[1][:300])

    def compare_naive(self, desc, t1, t2):
        self.total += 1
        f1, f2 = self.path('n1.nwk', t1), self.path('n2.nwk', t2)
        _, h = self.run(self.fable, f1, f2)
        _, nv = self.run(self.fable, '--naive', f1, f2)
        if len(h) != 2 or h != nv:
            self.bad += 1
            print('MISMATCH (vs naive)', desc, h, nv)

    def expect_same_rejection(self, desc, t1, t2):
        self.total += 1
        f1, f2 = self.path('e1.nwk', t1), self.path('e2.nwk', t2)
        frc, _ = self.run(self.fable, f1, f2)
        arc, _ = self.run(self.astra, f1, f2)
        if (frc != 0) != (arc != 0):
            self.bad += 1
            print('MISMATCH (rejection)', desc, frc, arc)


def random_cases(ck, seed, count, lo, hi):
    rng = random.Random(seed)
    shapes = ['bin', 'mix', 'high', 'cat', 'star']
    for i in range(count):
        n = rng.randint(lo, hi)
        if rng.random() < .1:
            labels = [rng.choice(['x%d', 'L_%d', "'s p%d'"]) % j for j in range(n)]
        else:
            labels = ['t%d' % j for j in range(n)]
        t1 = rtree(labels, rng.choice(shapes), rng)
        t2 = rtree(labels, rng.choice(shapes), rng)
        ck.compare('random seed=%d #%d n=%d' % (seed, i, n),
                   [decorate(ser(t, True), rng) + ';' for t in (t1, t2)],
                   [ser(t, False) + ';' for t in (t1, t2)])


def naive_cases(ck, seed, count):
    rng = random.Random(seed)
    for i in range(count):
        n = rng.randint(4, 14)
        labels = ['t%d' % j for j in range(n)]
        t1 = rtree(labels, rng.choice(['bin', 'mix', 'high']), rng)
        t2 = rtree(labels, rng.choice(['bin', 'mix', 'high']), rng)
        ck.compare_naive('naive seed=%d #%d n=%d' % (seed, i, n), ser(t1, True) + ';', ser(t2, True) + ';')


def big_cases(ck, n):
    rng = random.Random(5)
    a = ['t%d' % i for i in range(n)]
    b = a[:]
    rng.shuffle(b)
    def cat(l):
        c = l[0]
        for x in l[1:]:
            c = '(%s,%s)' % (c, x)
        return c + ';'
    cat_a, cat_b, star_a = cat(a), cat(b), '(' + ','.join(a) + ');'
    ck.compare('caterpillar vs caterpillar n=%d' % n, [cat_a, cat_b], [cat_a, cat_b])
    ck.compare('star vs caterpillar n=%d' % n, [star_a, cat_b], [star_a, cat_b])
    ck.compare('caterpillar vs star n=%d' % n, [cat_b, star_a], [cat_b, star_a])


def edge_cases(ck):
    same = [("a;", "a;"), ("(a,b);", "(a,b);"), ("(a,b);", "(b,a);"),
            ("((a,b),c);", "(a,(b,c));"), ("(a,b,c,d);", "((a,b),(c,d));"),
            ("((a,b),(c,d));", "((a,c),(b,d));")]
    for t1, t2 in same:
        ck.compare('edge %s | %s' % (t1, t2), [t1, t2], [t1, t2])
    for t1, t2 in [("(a,b,c);", "(a,b,d);"), ("(a,b,c);", "(a,b);"), ("(a,a,b);", "(a,a,b);")]:
        ck.expect_same_rejection('reject %s | %s' % (t1, t2), t1, t2)


def main():
    fable, astra = sys.argv[1], sys.argv[2]
    quick = '--quick' in sys.argv
    with tempfile.TemporaryDirectory() as tmp:
        ck = Checker(fable, astra, tmp)
        plan = [('edge cases', lambda: edge_cases(ck)),
                ('random n=4..12', lambda: random_cases(ck, 1, 100 if quick else 3000, 4, 12)),
                ('random n=4..40', lambda: random_cases(ck, 2, 100 if quick else 2000, 4, 40)),
                ('random n=4..300', lambda: random_cases(ck, 3, 30 if quick else 100, 4, 300)),
                ('random n=500..2500', lambda: random_cases(ck, 4, 2 if quick else 12, 500, 2500)),
                ('large caterpillar/star', lambda: big_cases(ck, 500 if quick else 3000)),
                ('vs naive oracle', lambda: naive_cases(ck, 99, 50 if quick else 300))]
        for name, fn in plan:
            before = ck.bad
            fn()
            print('%-26s %s' % (name, 'ok' if ck.bad == before else 'FAILED'), flush=True)
        print('%d cases, %d mismatches' % (ck.total, ck.bad))
        return 1 if ck.bad else 0


if __name__ == '__main__':
    result = []
    t = threading.Thread(target=lambda: result.append(main()))
    t.start()
    t.join()
    sys.exit(result[0] if result else 1)
