# Imports

import attitude_truth as attitude_truth


def main():
    
    print("Starting attitude determination...")
    
    attitude_truth.generate_truth_orbit()
    
    '''
    
    staticEstimation = calculateStaticAttitude()
    dynamicEstimation = calculateDynamicAttitude(staticEstimation)
    estimtedAttitude = runSynthesisFunction(staticEstimation, dynamicEstimation)
    
    return estimatedAttitude
    
    '''
    
    return 0


if __name__ == "__main__":
    main()
    
    