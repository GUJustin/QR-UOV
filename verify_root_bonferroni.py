#!/usr/bin/env python3
"""Numerical check of the explicit second-Bonferroni root lower bound."""
import math
Q=127
for r in (51,75,102):
    mu=1-1/Q
    pi=math.prod(1-Q**(-j) for j in range(1,r))
    bound=mu*pi/2-mu*mu/8+(Q-1)*mu*Q**(-r)/8
    derivative=pi/2-mu/4+(Q-1)*Q**(-r)/8
    assert derivative>0
    print(f"r={r}: p_root={bound:.15f}, derivative_at_mu0={derivative:.15f}")
