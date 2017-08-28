# -*- coding: utf-8 -*-

from wlcmcmc.wlcmcmclib import do_crank_shaft
import numpy as np
import pytest


@pytest.fixture
def straight_line_vertices():
    array = np.array([(0, 0, i) for i in range(6)], dtype=np.float64)
    yield array


class TestWLCMCMC:
    def test_do_crank_shaft(self, straight_line_vertices):
        actual = straight_line_vertices.copy()

        do_crank_shaft(actual, 1, 4, 2)
        expected = straight_line_vertices

        assert np.allclose(expected, actual)
