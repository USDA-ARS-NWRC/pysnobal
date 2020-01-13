#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pysnobal
----------------------------------

Tests for `pysnobal` module.
"""

import unittest

from pysnobal import pysnobal


class TestPysnobal(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_000_something(self):
        result = pysnobal.main()
        result


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
