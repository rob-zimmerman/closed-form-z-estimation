#!/usr/bin/env python3

"""Check large-shape gamma roots with 90-digit Decimal arithmetic.

The gamma simulation records the worst stable/direct root comparison in each
grid cell. This script recomputes the large-shape cases from the Bernoulli
series through order 30 and reports the largest absolute and relative errors.
"""

import csv
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 90

BERNOULLI = [
    (2, Decimal(1) / 6),
    (4, -Decimal(1) / 30),
    (6, Decimal(1) / 42),
    (8, -Decimal(1) / 30),
    (10, Decimal(5) / 66),
    (12, -Decimal(691) / 2730),
    (14, Decimal(7) / 6),
    (16, -Decimal(3617) / 510),
    (18, Decimal(43867) / 798),
    (20, -Decimal(174611) / 330),
    (22, Decimal(854513) / 138),
    (24, -Decimal(236364091) / 2730),
    (26, Decimal(8553103) / 6),
    (28, -Decimal(23749461029) / 870),
    (30, Decimal(8615841276005) / 14322),
]


def profile_score(shape):
    value = Decimal(1) / (2 * shape)
    for order, bernoulli in BERNOULLI:
        value += bernoulli / (Decimal(order) * shape**order)
    return value


def profile_score_prime(shape):
    value = -Decimal(1) / (2 * shape**2)
    for order, bernoulli in BERNOULLI:
        value -= bernoulli / shape ** (order + 1)
    return value


root = Path(__file__).resolve().parents[1]
audit_path = root / "output" / "summaries" / "gamma" / "numerical_accuracy.csv"
with audit_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

max_absolute = Decimal(0)
max_relative = Decimal(0)
checked = 0
for row in rows:
    target = Decimal(row["worst_profile_target"])
    stored = Decimal(row["stable_root"])
    report_as_dispersion = row["model"] == "Gamma regression"
    shape = Decimal(1) / stored if report_as_dispersion else stored
    if shape < 50:
        continue
    # Newton iteration is performed entirely in 90-digit Decimal arithmetic.
    for _ in range(12):
        shape -= (profile_score(shape) - target) / profile_score_prime(shape)
    reference = Decimal(1) / shape if report_as_dispersion else shape
    absolute = abs(stored - reference)
    relative = absolute / abs(reference)
    max_absolute = max(max_absolute, absolute)
    max_relative = max(max_relative, relative)
    checked += 1

print(f"checked={checked}")
print(f"max_absolute_difference={max_absolute:.8E}")
print(f"max_relative_difference={max_relative:.8E}")
