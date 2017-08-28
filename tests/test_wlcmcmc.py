# -*- coding: utf-8 -*-

import numpy as np
import pytest
import wlcmcmc.wlcmcmc as wlcmcmc
from wlcmcmc.vector import vnorm, edge_lengths
from collections import Counter


class WormLikeChainLib(wlcmcmc.WormLikeChain):
    USE_WLCMCMCLIB = True


class WormLikeChainPython(wlcmcmc.WormLikeChain):
    USE_WLCMCMCLIB = False


@pytest.fixture(scope='function')
def wlcs():
    wlc_lib = WormLikeChainLib(
        number_of_segments=100,
        segments_per_kuhn_length=10,
        kuhn_length=10.0,
        diameter=0.1,
        theta_max=1.3

    )
    wlc_python = WormLikeChainPython(
        number_of_segments=100,
        segments_per_kuhn_length=10,
        kuhn_length=10.0,
        diameter=0.1,
        theta_max=1.3

    )
    yield dict(wlc_lib=wlc_lib,
               wlc_python=wlc_python)


class TestWLCMCMC:
    def test_vertices_number(self, wlcs):
        for _, wlc in wlcs.items():
            assert 100 == len(wlc.vertices)

    def test_vertices_centered(self, wlcs):
        for _, wlc in wlcs.items():
            assert np.allclose(sum(wlc.vertices), [0.0])

    def test_vertices_distance_from_center(self, wlcs):
        for _, wlc in wlcs.items():
            assert np.allclose(vnorm(wlc.vertices), [wlc.edge_length / 2.0 / np.sin(2.0 * np.pi / 200)])

    def test_set_bending_rigidity(self, wlcs):
        # See Brownian Dynamics Simulation of Knot Diffusion along a Stretched DNA Molecule
        # Section DNA MODEL AND METHODS OF CALCULATIONS equation (2)
        # at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1367310/#__sec2title

        for _, wlc in wlcs.items():
            wlc.set_bending_rigidity(segments_per_kuhn_length=10)
            assert abs(4.81 - wlc.bending_rigidity) < 1e-2

    def test_edges(self, wlcs):
        for _, wlc in wlcs.items():
            # all lengths are one
            assert np.allclose(np.linalg.norm(wlc.edges, axis=1), [1])

            # edges point in all direction equally
            assert np.allclose(sum(wlc.edges), [0.0])

    def test_angles(self, wlcs):
        for _, wlc in wlcs.items():
            assert np.allclose(wlc.angles, [np.pi * 2.0 / 100])

    def test_get_energy(self, wlcs):
        for _, wlc in wlcs.items():
            e = wlc.get_energy()
            assert abs(0.5 * 4.81 * (np.pi * 2.0 / 100) ** 2 * 100 - e) < 1e-2

    def test_do_crankshaft(self, wlcs):
        def get_angle_counts(angles):
            return Counter([round(t, 4) for t in angles])

        for _, wlc in wlcs.items():
            assert Counter({0.0628: 100}) == get_angle_counts(wlc.angles)
            wlc.do_crankshaft(2, 10, np.pi / 2.0)
            assert Counter({0.0628: 98, 0.4429: 2}) == get_angle_counts(wlc.angles)

    def test_do_monte_carlo_steps(self, wlcs):
        for _, wlc in wlcs.items():
            angles = []

            # warm-up
            wlc.do_monte_carlo_steps(1000)

            for _ in range(100):
                wlc.do_monte_carlo_steps(100)
                angles.extend(list(wlc.angles))

            assert abs(0.55812197 - np.array(angles).mean()) < 2e-2
