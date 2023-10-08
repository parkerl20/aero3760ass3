# Imports

import attitude_truth as attitude_truth


def main():
    
    print("Starting attitude determination...")
    
    attitude_truth.generate_truth_orbit()
    
    return 0


if __name__ == "__main__":
    main()
    
    
    
    
    
'''

BRAINSTORM:

Static or dynamic attitude determination

Star Trackers onboard Satellite

Two vectors: taken from unit vectors to two stars tracked by two star trackers for fine attitude determination
We consider just the unit vectors, as the length of the vector has no information relevant to attitude determination

'''