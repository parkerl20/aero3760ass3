"""Rotation matrices for (3 x 1) vectors"""

import numpy as np


def rot3_x(theta: float, degrees: bool = False) -> np.ndarray:
    """Returns a 3x3 rotation matrix about the x-axis.
    
    Args:
        theta (float): The angle of rotation.
        degrees (bool, optional): Whether the angle is in degrees or radians. Defaults to False.
    
    Returns:
        np.ndarray: The X-rotation matrix.
    """
    theta = np.radians(theta) if degrees else theta
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)]
    ])

def rot3_y(theta: float, degrees: bool = False) -> np.ndarray:
    """Returns a 3x3 rotation matrix about the y-axis.

    Args:
        theta (float): The angle of rotation.
        degrees (bool, optional): Whether the angle is in degrees or radians. Defaults to False.

    Returns:
        np.ndarray: The Y-rotation matrix.
    """
    theta = np.radians(theta) if degrees else theta
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

def rot3_z(theta: float, degrees: bool = False) -> np.ndarray:
    """Returns a 3x3 rotation matrix about the z-axis.

    Args:
        theta (float): The angle of rotation.
        degrees (bool, optional): Whether the angle is in degrees or radians. Defaults to False.

    Returns:
        np.ndarray: The Z-rotation matrix.
    """
    theta = np.radians(theta) if degrees else theta
    
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])