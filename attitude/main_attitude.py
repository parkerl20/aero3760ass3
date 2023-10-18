from attitude import attitude_truth

def main():
    
    print("Starting attitude determination...")
    
    attitude_truth.generate_truth_orbit()
    
    '''
    
    estimtedAttitude = calculateStaticAttitude()
    
    return estimatedAttitude
    
    '''
    
    return 0


if __name__ == "__main__":
    main()