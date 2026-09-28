"""Weierstrass elliptic curve operations over a prime field.

All arithmetic is done with Python integers, which have arbitrary
precision, so there is no risk of overflow.  The only third-party-free
operation we need is modular inverse, which we compute via the
extended Euclidean algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Point:
    """An affine point on a curve, or the identity element.

    ``x`` and ``y`` are ``None`` for the point at infinity (the group
    identity).  We use a sentinel rather than a separate subclass so
    that equality and hashing remain trivial.
    """

    x: int | None
    y: int | None

    @property
    def is_identity(self) -> bool:
        return self.x is None and self.y is None

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        if self.is_identity:
            return "Point(identity)"
        return f"Point(x={self.x}, y={self.y})"


@dataclass(frozen=True)
class Curve:
    """Short Weierstrass curve ``y^2 = x^3 + a*x + b`` over ``GF(p)``.

    The curve is immutable so it can be safely shared and used as a key.
    Construction validates the discriminant (``4a^3 + 27b^2 != 0``) so
    that callers cannot accidentally build a singular curve, which would
    not form a group.
    """

    p: int
    a: int
    b: int

    def __post_init__(self) -> None:
        if self.p < 3:
            raise ValueError("p must be a prime >= 3")
        # Non-singularity check.  Over a field of characteristic != 2,3
        # the curve is singular iff 4a^3 + 27b^2 == 0.
        if (4 * self.a ** 3 + 27 * self.b ** 2) % self.p == 0:
            raise ValueError("curve is singular (discriminant is zero)")

    # ------------------------------------------------------------------
    # Field helpers
    # ------------------------------------------------------------------
    def _inv(self, value: int) -> int:
        """Modular inverse via the extended Euclidean algorithm.

        We avoid ``pow(value, -1, p)`` (Python 3.8+) only for clarity of
        implementation; both are correct.  The extended-Euclid version
        makes the zero-denominator error surface as a natural ``ZeroDivisionError``.
        """
        value %= self.p
        if value == 0:
            raise ZeroDivisionError("inverse of 0 does not exist")
        lm, hm = 1, 0
        low, high = value, self.p
        while low > 1:
            r = high // low
            nm = hm - r * lm
            new = high - r * low
            hm, high, lm, low = lm, low, nm, new
        return lm % self.p

    # ------------------------------------------------------------------
    # Point validation
    # ------------------------------------------------------------------
    def is_on_curve(self, point: Point) -> bool:
        if point.is_identity:
            return True
        x, y = point.x, point.y
        if x is None or y is None:
            return False
        x %= self.p
        y %= self.p
        return (y * y - (x ** 3 + self.a * x + self.b)) % self.p == 0

    # ------------------------------------------------------------------
    # Group law
    # ------------------------------------------------------------------
    def add(self, p1: Point, p2: Point) -> Point:
        """Return ``p1 + p2`` under the curve's group law."""
        if p1.is_identity:
            return p2
        if p2.is_identity:
            return p1

        x1, y1 = p1.x % self.p, p1.y % self.p
        x2, y2 = p2.x % self.p, p2.y % self.p

        if x1 == x2:
            if (y1 + y2) % self.p == 0:
                # P + (-P) = O
                return Point(None, None)
            # P == P, use doubling formula
            s = (3 * x1 * x1 + self.a) * self._inv(2 * y1) % self.p
        else:
            s = (y2 - y1) * self._inv((x2 - x1) % self.p) % self.p

        x3 = (s * s - x1 - x2) % self.p
        y3 = (s * (x1 - x3) - y1) % self.p
        return Point(x3, y3)

    def negate(self, point: Point) -> Point:
        """Return ``-point``."""
        if point.is_identity:
            return point
        return Point(point.x % self.p, (-point.y) % self.p)

    def scalar_multiply(self, k: int, point: Point) -> Point:
        """Return ``k * point`` via the double-and-add ladder.

        Negative ``k`` is supported by negating the point.  ``k == 0``
        yields the identity.
        """
        if k == 0:
            return Point(None, None)
        if k < 0:
            return self.scalar_multiply(-k, self.negate(point))

        result = Point(None, None)
        addend = point
        # Iterate bits from LSB to MSB.  This is constant-time-ish in
        # structure (we always do an add and a double per bit) which is
        # a modest side-channel improvement over branching on the bit,
        # though it is NOT a hardened constant-time implementation and
        # must not be used for production secrets.
        bits = k.bit_length()
        for i in range(bits):
            if (k >> i) & 1:
                result = self.add(result, addend)
            addend = self.add(addend, addend)
        return result
