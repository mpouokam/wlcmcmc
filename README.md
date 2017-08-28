# Worm-Like Chain, Markov Chain Monte Carlo

Simulate closed Worm-Like Chains (Kratky-Porod model) like no other!

To install, activate the appropriate virtual environment and then:

```make install```

`make install` runs SWIG, installs NumPy if necessary, and then installs this package.

## Generating a Worm-Like Chain

Generating WLC conformations is easy. For instance if you wanted to find the angle distribution for a certain set of
paramters you could do the following:

```
wlc = WormLikeChainLib(
    number_of_segments=100,
    segments_per_kuhn_length=10,
    kuhn_length=10.0,
    diameter=0.1,
    theta_max=1.3
)

# warm-up
wlc.do_monte_carlo_steps(1000)

angles = []

for _ in range(100):
    wlc.do_monte_carlo_steps(100)
    angles.extend(list(wlc.angles))
```