import unittest

from elliptic_curve_operations import Curve, Point


# secp256k1 parameters (the Bitcoin curve) — a real, well-known curve
# that gives us non-trivial field size to exercise.
SECP256K1_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP256K1_A = 0
SECP256K1_B = 7

# A generator point on secp256k1.
G = Point(
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)


class TestCurveConstruction(unittest.TestCase):
    def test_rejects_singular_curve(self):
        # y^2 = x^3 over GF(5): 4*0 + 27*0 == 0
        with self.assertRaises(ValueError):
            Curve(5, 0, 0)

    def test_rejects_small_prime(self):
        with self.assertRaises(ValueError):
            Curve(2, 1, 1)

    def test_accepts_valid_curve(self):
        # y^2 = x^3 + x + 1 over GF(5): 4*1 + 27*1 = 31 = 1 (mod 5) != 0
        c = Curve(5, 1, 1)
        self.assertEqual(c.p, 5)


class TestPoint(unittest.TestCase):
    def test_identity_property(self):
        self.assertTrue(Point(None, None).is_identity)
        self.assertFalse(Point(1, 2).is_identity)


class TestIsOnCurve(unittest.TestCase):
    def setUp(self):
        self.curve = Curve(SECP256K1_P, SECP256K1_A, SECP256K1_B)

    def test_generator_is_on_curve(self):
        self.assertTrue(self.curve.is_on_curve(G))

    def test_identity_is_on_curve(self):
        self.assertTrue(self.curve.is_on_curve(Point(None, None)))

    def test_random_point_not_on_curve(self):
        self.assertFalse(self.curve.is_on_curve(Point(1, 1)))

    def test_point_with_none_components_not_on_curve(self):
        # Only (None, None) is the identity; (1, None) is malformed.
        self.assertFalse(self.curve.is_on_curve(Point(1, None)))


class TestAddition(unittest.TestCase):
    def setUp(self):
        self.curve = Curve(SECP256K1_P, SECP256K1_A, SECP256K1_B)

    def test_identity_is_additive_identity(self):
        O = Point(None, None)
        self.assertEqual(self.curve.add(G, O), G)
        self.assertEqual(self.curve.add(O, G), G)
        self.assertEqual(self.curve.add(O, O), O)

    def test_point_plus_its_negation_is_identity(self):
        neg_g = self.curve.negate(G)
        self.assertEqual(self.curve.add(G, neg_g), Point(None, None))

    def test_doubling(self):
        two_g = self.curve.add(G, G)
        expected = Point(
            0xC6047F9441ED7D6D3045406E95C07CD85C778E4B8CEF3CA7ABAC09B95C709EE5,
            0x1AE168FEA63DC339A3C58419466CEAEEF7F632653266D0E1236431A950CFE52A,
        )
        self.assertEqual(two_g, expected)
        self.assertTrue(self.curve.is_on_curve(two_g))

    def test_addition_of_distinct_points(self):
        two_g = self.curve.add(G, G)
        three_g = self.curve.add(two_g, G)
        # 3G on secp256k1
        expected = Point(
            0xF9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9,
            0x388F7B0F632DE8140FE337E62A37F3566500A99934C2231B6CB9FD7584B8E672,
        )
        self.assertEqual(three_g, expected)

    def test_negate_identity(self):
        self.assertEqual(self.curve.negate(Point(None, None)), Point(None, None))


class TestScalarMultiply(unittest.TestCase):
    def setUp(self):
        self.curve = Curve(SECP256K1_P, SECP256K1_A, SECP256K1_B)

    def test_zero_times_point_is_identity(self):
        self.assertEqual(self.curve.scalar_multiply(0, G), Point(None, None))

    def test_one_times_point_is_point(self):
        self.assertEqual(self.curve.scalar_multiply(1, G), G)

    def test_two_via_scalar_matches_addition(self):
        self.assertEqual(
            self.curve.scalar_multiply(2, G),
            self.curve.add(G, G),
        )

    def test_known_nG_value(self):
        # n is the order of the secp256k1 group; n*G == identity.
        n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
        self.assertEqual(self.curve.scalar_multiply(n, G), Point(None, None))

    def test_negative_scalar(self):
        # k*G and (-k)*G should be negatives of each other.
        pos = self.curve.scalar_multiply(7, G)
        neg = self.curve.scalar_multiply(-7, G)
        self.assertEqual(self.curve.add(pos, neg), Point(None, None))

    def test_distributivity(self):
        # (a+b)*G == a*G + b*G
        a, b = 11, 17
        lhs = self.curve.scalar_multiply(a + b, G)
        rhs = self.curve.add(
            self.curve.scalar_multiply(a, G),
            self.curve.scalar_multiply(b, G),
        )
        self.assertEqual(lhs, rhs)


class TestSmallCurve(unittest.TestCase):
    """Exhaustive-ish tests on a tiny curve where we can reason by hand."""

    def setUp(self):
        # y^2 = x^3 + 2x + 2 over GF(17).  This is non-singular.
        self.c = Curve(17, 2, 2)
        # (5, 1) is on this curve: 1 == 125 + 10 + 2 = 137 == 1 (mod 17)
        self.P = Point(5, 1)

    def test_point_on_small_curve(self):
        self.assertTrue(self.c.is_on_curve(self.P))

    def test_double_on_small_curve(self):
        # 2P computed by hand: slope = (3*25 + 2) / (2*1) = 77 / 2 = 9 * inv(2)
        # inv(2, 17) = 9, so s = 9*9 = 81 = 13 (mod 17)
        # x3 = 13^2 - 2*5 = 169 - 10 = 159 = 6 (mod 17)
        # y3 = 13*(5 - 6) - 1 = -13 - 1 = -14 = 3 (mod 17)
        two_p = self.c.add(self.P, self.P)
        self.assertEqual(two_p, Point(6, 3))
        self.assertTrue(self.c.is_on_curve(two_p))

    def test_scalar_consistent_with_repeated_add(self):
        acc = Point(None, None)
        for _ in range(19):
            acc = self.c.add(acc, self.P)
        self.assertEqual(self.c.scalar_multiply(19, self.P), acc)


if __name__ == "__main__":
    unittest.main()
