import unittest

from koffi import PotentialSource

class TestImageMetadata(unittest.TestCase):
    def test_init(self):
        ps = PotentialSource()

        self.assertEqual(ps.position_at, {})
        self.assertEqual(ps.times, None)

    def test_build_from_times_and_known_positions(self):
        pos = [
            [200.501433, -14.166194],
            [200.502433, -14.166194],
            [200.503433, -14.166194],
            [200.504433, -14.166194],
            [200.505433, -14.166194]
        ]

        mjds = [
            57131.25239706122,
            57131.25888067302,
            57133.24996839132,
            57133.25643622691,
            57133.26293793291
        ]

        ps = PotentialSource()
        ps.build_from_times_and_known_positions(pos, mjds)

        self.assertEqual(ps[57133.25643622691][0], 200.504433)

if __name__ == '__main__':
    unittest.main()