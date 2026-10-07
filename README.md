# elliptic-curve-operations

Point addition and scalar multiplication on short Weierstrass curves `y^2 = x^3 + a*x + b` over a prime field `GF(p)`, using only the Python standard library.

## Usage

```python
from elliptic_curve_operations import Curve, Point

# secp256k1 (the Bitcoin curve)
curve = Curve(
    p=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F,
    a=0,
    b=7,
)

G = Point(
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)

assert curve.is_on_curve(G)

two_g = curve.add(G, G)
secret = 0x18E14A7B6A307F426A94F8114701E7C8E774E7F9A47E2C2035DB29A206321725
pub = curve.scalar_multiply(secret, G)
print(pub)
```

Exported names: `Curve` (frozen dataclass with fields `p`, `a`, `b`), `Point` (frozen dataclass with fields `x`, `y`, either `int` or `None` for the identity).

`Curve` methods: `is_on_curve(point)`, `add(p1, p2)`, `negate(point)`, `scalar_multiply(k, point)`.

## Why this exists

For cryptographic prototyping and teaching, you often need elliptic curve group operations without pulling in a heavy dependency. This library implements the affine-coordinate group law directly, with no C extensions and no third-party packages. The trade-off: affine coordinates are slower than projective/Jacobian for scalar multiplication, and the implementation is not hardened against side-channel attacks — it is for prototyping and education, not for guarding real secrets.

## Awkward edges

- The point at infinity (group identity) is represented as `Point(None, None)`. A `Point` with only one coordinate `None` is considered malformed and `is_on_curve` returns `False` for it.
- `Curve` rejects singular curves (where `4a^3 + 27b^2 ≡ 0 mod p`) at construction time, since those do not form a group.
- `scalar_multiply` accepts negative `k` by negating the point; `k=0` returns the identity.
- The double-and-add ladder always performs both an add and a double per bit, which is a mild side-channel improvement over bit-branching, but this is still not a constant-time implementation. Do not use it for production keys.

## Design notes

The window stores values eagerly rather than keeping running aggregates. Running
sums drift with floating point over long streams, and recomputing from a small
buffer is cheap enough that the drift is not worth the speed.

## Limitations

Values are coerced to floats, so very large integers lose precision. If you need
exact integer aggregates over a window, this is the wrong tool.

