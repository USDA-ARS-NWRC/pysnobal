#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pysnobal
----------------------------------

Tests for `pysnobal` module.
"""

import unittest

from pysnobal.snobal import PySnobal


class TestPysnobal(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_000_something(self):
        result = PySnobal().run()
        result
