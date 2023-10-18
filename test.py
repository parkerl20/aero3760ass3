from spacesim import satellite as sat
from spacesim import orbit as orb
from spacesim import celestial_body as cb
from spacesim import sensor
from spacesim import estimation as est

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime as dt
import numpy as np

def imu_simulator(
    imu: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    imu.mesurement = np.array([1, 2, 3])
    return True

def ekf_simulator(
    time: float,
    ekf: sat.SatelliteAlgorithm,
    satellite: sat.RealTimeSatellite,
    sensors: dict[str, sensor.SatelliteSensor]
) -> None:
    imu = sensors.get("imu", None)
    
    if imu is not None:
        imu_measurement = imu.get_measurement()
        print(f"EKF: {imu_measurement} at {time}")
    else:
        print(f"EKF: No measurement at {time}")

def main() -> None:
    a = 6932386.57371916
    e = 0
    i = -33
    arg_p = 0
    RAAN = 238.82
    theta = 0
    epoch = dt.datetime(2023, 1, 1)
    
    propagation_time = 20
    
    SPY = orb.Orbit(
        a,
        e,
        i,
        RAAN,
        arg_p,
        theta,
        cb.earth,
        epoch,
        name="SPY",
        colour="blue"
    )
    
    realtime_sat = sat.RealTimeSatellite(
        a,
        e,
        i,
        RAAN,
        arg_p,
        theta,
        cb.earth,
        epoch,
        name="Real Time",
        colour="red",
        propagation_length=10,
        propagation_step=1
    )
    
    
    imu_sensor = sensor.SatelliteSensor(
        "imu",
        np.eye(3),
        imu_simulator,
        frequency=0.5
    )
    
    ekf = est.ExtendedKalmanFilter(
        None,
        np.zeros(6),
        np.eye(6)
    )
    
    ekf_algorithm = sat.SatelliteAlgorithm(
        "Orbit Determination EKF",
        ekf,
        ekf_simulator
    )
    
    realtime_sat.attach_sensor(imu_sensor)
    realtime_sat.add_algorithm(ekf_algorithm)
    
    # r_spy, v_spy, t_spy = SPY.propagate(
    #     propagation_time,
    #     max_step=1
    # )
    
    # DO real time orbit
    r_rt, v_rt, t_rt = [], [], []
    
    for r_i, v_i, t_i in realtime_sat:
        if t_i > propagation_time:
            break
    
        print(t_i)
        r_rt.append(r_i)
        v_rt.append(v_i)
        t_rt.append(t_i)
        # print(f"Real time: {t_i}")
        
    
    r_rt = np.array(r_rt).T
    v_rt = np.array(v_rt).T
    t_rt = np.array(t_rt)
    
    # Plot
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    
    # ax.plot(r_spy[0], r_spy[1], r_spy[2], label="SPY")
    # ax.plot(r_rt[0], r_rt[1], r_rt[2], label="Real Time")
    # ax.legend()
    
    # plt.show()
    

if __name__ == "__main__":
    main()