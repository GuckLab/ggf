import numpy as np


def get_stress(sigma_0, order=2, n_points=100):
    """sigma_0 cos^n(theta) stress model"""
    theta = np.linspace(0, np.pi, n_points, endpoint=True)
    sigma = sigma_0 * (np.cos(theta))**order
    return sigma
