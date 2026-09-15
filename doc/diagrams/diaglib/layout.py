"""Tiny layout helpers for the box-and-arrow figures."""


def spread(x0, x1, n, w):
    """Left edges of ``n`` boxes of width ``w`` spread evenly over ``[x0, x1]``."""
    if n == 1:
        return [x0 + (x1 - x0 - w) / 2]
    gap = (x1 - x0 - n * w) / (n - 1)
    return [x0 + i * (w + gap) for i in range(n)]


def stack(y0, n, h, gap):
    """Top edges of ``n`` boxes of height ``h`` stacked from ``y0``."""
    return [y0 + i * (h + gap) for i in range(n)]
