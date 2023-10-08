# Imports

import attitude_truth as attitude_truth


def main():
    
    print("Starting attitude determination...")
    
    attitude_truth.generate_truth_orbit()
    
    '''
    
    calculateStaticAttitude()
    calculateDynamicAttitude()
    estimtedAttitude = runSynthesisFunction()
    
    '''
    
    return 0


if __name__ == "__main__":
    main()
    
    