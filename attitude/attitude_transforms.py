import numpy as np



def euler2quat(x):
    """
    Converts Euler angles to quaternion
    Inputs:
        x: Euler angles in degrees yaw, pitch, roll
            in order (z,y,x: yaw, pitch, roll: psi, theta, phi)
    Outputs:
        q: normalised quaternion in form w x y z
    """

    yaw = np.deg2rad(x[0])
    pitch = np.deg2rad(x[1])
    roll = np.deg2rad(x[2])
    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    
    # normalise
    quat = np.array([ qw,qx, qy, qz])
    quat = quat/np.linalg.norm(quat)
    return quat

def quat2euler(q):
    """
    Converts quaternions to Euler angles
    Inputs:
        q: normalised quaternion in form w x y z
        
    Outputs:
        x: Euler angles in degrees yaw, pitch roll
            in order (z,y,x: yaw, pitch, roll: psi, theta, phi)
        
    """
    # roll (x-axis rotation)
    #x,y,z,w
    (w, x, y, z) = (q[0], q[1], q[2], q[3])

    t0 = 2 * (w * x + y * z)
    t1 = 1 - 2 * (x * x + y * y)
    roll = np.arctan2(t0, t1)
    t2 = 2 * (w * y - z * x)
    t2 = 1 if t2 > 1 else t2
    t2 = -1 if t2 < -1 else t2
    pitch = np.arcsin(t2)
    t3 = 2 * (w * z + x * y)
    t4 = 1 - 2 * (y * y + z * z)
    yaw = np.arctan2(t3, t4)

    yaw_deg = np.rad2deg(yaw)
    pitch_deg = np.rad2deg(pitch)
    roll_deg = np.rad2deg(roll)
    return np.array([yaw_deg, pitch_deg, roll_deg])

def quaternion_multiply(q1, q2):
    """
    Multiplies 2  quaternions (do not need to be normalised!):
    performs q2 rotation first then q1

    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2

    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2

    return np.array([w, x, y, z])

